import sys
import os
import re
import itertools
import argparse
import logging
import numpy as np
import numba
import pickle
import h5py
import toml
import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, '../utils')
sys.path.insert(0, utils_path)

import hdf_utils

# NOTE: adjust to produce a dictionary, keys of which are chromosome
#   identifiers, values are current numpy arrays containing coverage
#   and bootstrapped coverage.

def progress(count, total, status=''):
    '''Print progress to stdout. 
    Taken from https://gist.github.com/vladignatyev/06860ec2040cb497f0f3
    '''
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = count / float(total)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[{}] {:.1%} ...{}\n'.format(bar, percents, status))
    sys.stdout.flush()

def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks. Taken from
    itertools documentation

    Test:
    >>> out = grouper('ABCDEFG', 3, 'x')
    >>> for group in out:
    ...    print group
    ('A', 'B', 'C')
    ('D', 'E', 'F')
    ('G', 'x', 'x')
    """
    args = [iter(iterable)] * n
    return itertools.zip_longest(fillvalue=fillvalue, *args)

class SamAlignment:
    """ Class for holding a single line out of a SAM file. Along with methods
    to manipulate and extract data from the alignment.
    """
        
    def __init__(self, line, sep = '\t'):
        """ Class initalizer. Takes in a single line from a Sam alignment and
        stores each and every field as a class attribute.

        Input: line - str Raw Sam input file line.
        Modifies: self.QNAME - str Name of the read
                  self.FLAG  - int Bitwise flag of read information
                  self.RNAME - str Name of the reference chromosome
                  self.POS   - int 0-based left most mapping of first base
                  self.MAPQ  - int mapping quality
                  self.CIGAR - str contains information on mapping with gaps
                  self.RNEXT - str name of reference chromosome of read pair
                  self.PNEXT - int POS of other read in pair
                  self.TLEN  - int number of bases from left to right in pair
                  self.SEQ   - str contains sequence of the read
                  self.QUAL  - str contains ASCII base quality information
                  self.OPT   - str optional fields specific to each aligner
        Returns: None

        Tests:

        Test reading in a typical SAM entry
        >>> line = "testQ.1.testR.20.30.3M.testR.25.9.AAA.(((.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.QNAME == "testQ"
        True
        >>> read.FLAG == 1
        True
        >>> np.all(read.BITS == [1,0,0,0,0,0,0,0,0,0,0,0])
        True
        >>> read.RNAME == "testR"
        True
        >>> read.POS == 19
        True
        >>> read.MAPQ == 30
        True
        >>> read.CIGAR == "3M"
        True
        >>> read.RNEXT == "testR"
        True
        >>> read.PNEXT == 24
        True
        >>> read.TLEN == 9
        True
        >>> read.SEQ == "AAA"
        True
        >>> read.QUAL == "((("
        True
        >>> read.OPT["NM"] == 0
        True
        """

        linearr = line.rstrip().split(sep)
        # All information comes from SAM Specificationv1 18 Nov 2015
        # This is the name of the read, two reads from the same
        # template share the same QNAME. * indicates the information
        # is unknown.
        self.QNAME = linearr[0]
        # This is a combination of bitwise flags. Each bit is as follows:
        # Bit    Hex    Description
        # 1      0x1    template having multiple segments
        # 2      0x2    each segment properly aligned
        # 4      0x4    segment unmapped
        # 8      0x8    next segment in template unmapped
        # 16     0x10   SEQ being reverse complemented
        # 32     0x20   SEQ of the next segment reverse complemented
        # 64     0x40   the first segment in the template
        # 128    0x80   the last segment in the template
        # 256    0x100  secondary alignment
        # 512    0x200  not passing filters(platform/vendor)
        # 1024   0x400  PCR or optical duplicate
        # 2048   0x800  supplementary alignment
        self.FLAG = int(linearr[1])
        self.BITS = np.zeros(12, dtype=uint8)
        bits = f"{self.FLAG:b}"
        bit_arr = [ int(_) for _ in bits[::-1] ]
        for i,bit in enumerate(bit_arr):
            self.BITS[i] = bit
        # Reference sequence name. Would be name of the chromosome in an
        # organism with multiple chromosomes
        self.RNAME = linearr[2]
        # 0-based left most mapping of the first base. -1 if unmapped. If
        # -1 no assumptions can be made about RNAME and CIGAR
        self.POS = int(linearr[3]) - 1
        # Mapping quality. Equal to -10log10(P(mapping position is wrong)).
        # Rounded to nearest integer. 255 means it is unavailable
        self.MAPQ = int(linearr[4])
        # CIGAR string, set to * if unavailable.
        # Value  Description
        #   M     alignment match
        #   I     insertion to the reference
        #   D     deletion from the reference
        #   N     skipped region from the reference
        #   S     soft clipping (clipped sequences present in SEQ)
        #   H     hard clipping (clipped sequences NOT in SEQ)
        #   P     padding (silent deletion from padded reference)
        #   =     sequence match
        #   X     sequence mismatch
        # H can only be present as first or last operation
        # S may only have H operations between them and the ends
        # N operation represents an intron for mRNA-Genome
        # Sum of the lengths shall equal the length of SEQ
        self.CIGAR = linearr[5]
        # Reference sequence name for other read in template, set to
        # = if they are identical and * if unavailable.
        self.RNEXT = linearr[6]
        # POS of the other read in template. -1 when information is unavail.
        # If -1 no assumptions can be made on RNEXT and bit 0x20
        # Subtract one to keep the value in 0-based coordinates
        self.PNEXT = int(linearr[7]) - 1
        # If all mapped to same reference. Then equal number of bases
        # from leftmost mapped base to righmost mapped base. Set to 0
        # for single-segment template or when information is unavail.
        self.TLEN = int(linearr[8])
        # Sequence for the segment
        self.SEQ = linearr[9]
        # ASCII quality score plus 33 (Phred 33)
        self.QUAL = linearr[10]
        # Optional fields. See SAM spec for more details
        self.OPT = self.parse_opt_fields(linearr[11:])

        ## Attributes that act as caches
        self.cigar_tuples = None
        self.aligned_blocks = None
        self.aligned_seq = None
        self.aligned_phred = None
        self.aligned_locs = None
        self.aligned_reference = None
        self.aligned_muts = None

    def parse_opt_fields(self, fields_list):
        """ Parses each field from the optional fields in the 
        SAM file. 

        Inputs: List. Contains all the optional fields in a sam file
        Modifies: nothing
        Returns: dictionary. Has field tags as keys and correctly converted
                 values. 'H' and 'B' field types are converted to strings
                 for now. Could be implemented to arrays in the future
        Tests:
        >>> line = "testQ.4.testR.20.30.3M.testR.25.9.AAA.(((.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> ans = read.parse_opt_fields(["NM:i:0", "JK:f:0", "HT:Z:0", "QR:A:0",
        ...                              "AB:H:0", "CD:B:0"])
        >>> type(ans["NM"]) is int
        True
        >>> type(ans["JK"]) is float
        True
        >>> type(ans["HT"]) is str
        True
        >>> type(ans["QR"]) is str
        True
        >>> type(ans["AB"]) is str
        True
        >>> type(ans["CD"]) is str
        True
        """
        # Technically, H and B are arrays, but I don't have that implemented
        # right now, could be implemented in the future
        # Have a dictionary to hold the functions to convert each type of
        # field to the type it indicates
        d_type_func = {'A': str, 'i': int, 'f':float, 'Z':str, 
                       'H': str, 'B': str}
        # initialize a dictionary to hold the parsed fields
        field_dict = {}
        # loop through the additional fields
        for field in fields_list:
            # get name, type and value
            name, d_type, value  = field.split(':', 3)
            # convert value to appropriate type and put in dictionary for
            # holding fields
            field_dict[name] = d_type_func[d_type](value)
        # return this dictionary
        return(field_dict)

    def get_cigar_tuples(self):

        """Function to parse the cigar field and turn it into tuples
        using regular expressions. Will also cache the return value into
        self.cigar_tuples so that the calculation only has to be performed
        once, when handling the record.

        Input: nothing
        Modifies: self.cigar_tuples
        Returns: a tuple of the cigar field broken up by separators.

        Tests:

        All M test
        >>> line = "testQ.1.testR.20.30.7M.testR.25.9.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_cigar_tuples()
        [(7, 'M')]

        More complex test
        >>> read.CIGAR = "3M4D2M1I1M" 
        >>> read.cigar_tuples = None
        >>> read.get_cigar_tuples()
        [(3, 'M'), (4, 'D'), (2, 'M'), (1, 'I'), (1, 'M')]

        """
        # if the value returned by this function has already been cached, then
        # just return that value and skip the calculation
        if self.cigar_tuples:
            return self.cigar_tuples

        # otherwise initalize a list to store each of the tuples
        cigar_tuples = []
        # split using regular expressions. The parenthesis in the re gives us
        # the seperators as well
        cigar_list = re.split('([MIDNSHP=X])', self.CIGAR)
        
        # loop through by twos using itertools grouper recipe from the 
        # python itertools documentation.
        # The cigar string always starts with a number and ends in a char,
        # so we cut the list short by one since the last value by twos will
        # end up being a None.
        for number, char in grouper(cigar_list[:-1], 2):
            cigar_tuples.append((int(number), char))
        # set the value into the cache so the calculation doesn't have to be
        # performed again
        self.cigar_tuples = cigar_tuples
        # return the now cached cigar_tuples value
        return self.cigar_tuples

    def get_aligned_blocks(self):

        """ Function to take the cigar field and determine the locations
        where there is continuous mapping coverage. 
        
        Returns: a list of [start, end, strand] locations where continuous
                 mapping coverage occurs. Continuous coverage includes
                 locations where there is continuous 'M', '=', 'D', or 'X'
                 in the CIGAR field. 
                 Breaks occur at 'S', 'I', or 'N' in the CIGAR field.

        Tests:
        All matching test, read.POS is at 1, cigar is 7M
        >>> line = "testQ.1.testR.2.30.7M.testR.4.9.AGTCGCT.!#%()+,.NM:i:0"
        >>> read = SamAlignment(line, sep = '.')
        >>> read.get_aligned_blocks()
        [(1, 8)]

        Test internal deletions
        >>> read.CIGAR = "2M3D5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_aligned_blocks()
        [(1, 11)]

        Test internal introns
        >>> read.CIGAR = "2M3N5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_aligned_blocks()
        [(1, 3), (6, 11)]

        Test internal insertions
        >>> read.CIGAR = "2M3I2M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_aligned_blocks()
        [(1, 5)]

        Test soft clipping on either end
        >>> read.CIGAR = "2S5M"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_aligned_blocks()
        [(1, 6)]
        >>> read.CIGAR = "5M2S"
        >>> read.cigar_tuples = None
        >>> read.aligned_blocks = None
        >>> read.get_aligned_blocks()
        [(1, 6)]
        """
        if self.aligned_blocks:
            return self.aligned_blocks
        # initialize variables needed in the loop
        end = None
        pos_list = []
        if self.BITS[4]:
            strand = 1
        else:
            strand = 0
        start = self.POS
        # parse the cigar string into tuples
        cigar_tuples = self.get_cigar_tuples()
       
        # go through each cigar number and character
        for val, char in cigar_tuples:
            # if it is an alignment match of any sort,
            # just add the value to determine the end. If the end
            # had not been determined previously, add it. otherwise,
            # just add the value to the previous end.
            if char in ['M', '=', 'X', 'D']:
                if end is not None:
                    end += val 
                else:
                    end = start + val
            # If there is an intron, go ahead and append
            # the previous [start, end, strand] list to the list and reset the
            # next end to being unknown. Additionally, push the next start
            # to the end of the intron. If an end had not been 
            # determined yet, do not add the (start, end) pair to the list
            # as this is a continuing insertion of some sort.
            elif char == 'N':
                if end is not None:
                    pos_list.append([start, end, strand])
                    start = end + val
                    end = None
                else:
                    start += val
            elif char in ['S', 'I']:
                continue
        # Finally, once all the way through the cigar field. Append the last
        # [start, stop, strand] pair as long as it doesn't end in an insertion. If
        # it ends in an insertion then don't append the last pair.
        if end is not None:
            pos_list.append([start, end, strand])

        # last add the final list to the cache in the class and return
        # the value.

        self.aligned_blocks = pos_list
        return self.aligned_blocks


class ReadSampler(object):
    def __init__(self):
        self.reads = []
        self.total = 0
        self.sampling = False

    def add_read(self, new_read):
        '''Adds a "new_read" to this ReadSampler object.
        
        Args:
        -----
        new_read : tuple
            new_read contains the following: (start,end,strand,rname),
            where start is the read's start position, end is the read's
            end position, and rname is the reference identifier to which
            the read aligned. rname would be the reference's unique identifier
            for genomes with multiple dna elements.
            For the strand element, it is 0 if R1 read aligned to "+" strand,
            and 1 if R1 read aligned to "-" strand.

        Modifies:
        ---------
        self.reads : list
            appends new_read to the reads list.
        self.total : int
            adds one to the total number of reads.
        '''
        if self.sampling:
            self.convert_to_list()
        self.reads.append(new_read)
        self.total += 1

    def add_reads(self, new_reads):
        if self.sampling:
            self.convert_to_list()
        self.reads.extend(new_reads)
        self.total += len(new_reads)
    
    def convert_to_array(self):
        self.reads = np.asarray(self.reads, dtype=np.uint64)
        self.sampling = True

    def convert_to_list(self):
        self.reads = list(self.reads)
        self.sampling = False

    def pull_read(self):
        if not self.sampling:
            self.convert_to_array()
        index = np.random.randint(0, self.total)
        return self.reads[index, :]

    def pull_reads(self, n):
        if not self.sampling:
            self.convert_to_array()
        index = np.random.randint(0, self.total, size=n)
        index = np.sort(index)
        return self.reads[index, :]

    def sort_reads(self):
        if not self.sampling:
            self.convert_to_array()
        # sort reads by chrom and pos
        sort_inds = np.lexsort((
            self.reads[:,0],
            self.reads[:,2]
        ))
        self.reads = self.reads[sort_inds,:]

    def from_hdf5(self, hdf_name):
        try:
            with h5py.File(hdf_name, 'r') as hf:
                self.reads = hf["parser"][...]
        except KeyError as e:
            print("==============================================")
            print(
                f"The hdf5 file {hdf_name} has no 'parser' dataset. This "\
                f"usually indicates that either read preprocessing or "\
                f"alignment failed. Check your fastq and bam files, "\
                f"and check your cutadapt, trimmomatic, and bowtie log "\
                f"and err files."
            )
            print("==============================================")
        self.total = self.reads.shape[0]
        self.sampling = True

def merge(intervals, strand):
    intervals.sort(key=lambda x: x[0])
    # take the first interval
    merged = [intervals[0]]
    # loop through all the intervals
    for this_interval in intervals:
        # if this interval starts within current merged interval, do the following
        if this_interval[0] <= merged[-1][1]:
            merged[-1] = [
                merged[-1][0], # same start
                max(merged[-1][1], this_interval[1]), # max end
                strand, # set strand
            ]
        else:
            merged.append(this_interval)
    return merged

def get_paired_blocks(r1, r2):
    if r1.BITS[4]:
        strand = 1
    else:
        strand = 0
    if r1.TLEN > 0 and r2.TLEN < 0:
        # each of left, right are (start,end,strand) lists
        left = r1.get_aligned_blocks()
        right = r2.get_aligned_blocks()
    elif r1.TLEN < 0 and r2.TLEN > 0:
        left = r2.get_aligned_blocks()
        right = r1.get_aligned_blocks()
    elif r1.POS == r2.POS:
        left = r1.get_aligned_blocks()
        right = r2.get_aligned_blocks()
    else:
        raise ValueError("Pair not consistent {} {}".format(r1.QNAME, r2.QNAME))
    total_blocks = []
    # if the right-most position in "left"
    # is less than the left-most position in "right"
    # do the following
    if left[-1][1] < right[0][0]:
        total_blocks.append([left[-1][1], right[0][0]])
    total_blocks.extend(left)
    total_blocks.extend(right)
    total_blocks = merge(total_blocks, strand)
    if len(total_blocks) > 1:
        raise RuntimeError("Gapped read found {} {} {}".format(
            r1.QNAME, r2.QNAME, str(total_blocks))
        )
    return total_blocks[0]

def create_read_list(samfile, ctg_lut):
    read_sampler = ReadSampler()
    for line in samfile:
        line = SamAlignment(line)
        # vals is [[start, end, strand]]
        vals = line.get_aligned_blocks()
        ctg_idx = ctg_lut[line.RNAME]["idx"]
        if len(vals) > 1:
            logging.info("Skipping gapped read {} {}".format(
                line.QNAME, str(vals)
            ))
        vals = vals[0]
        vals.append(ctg_idx)
        read_sampler.add_read(vals)
    return read_sampler

def create_read_list_paired(samfile, ctg_lut):
    read_sampler = ReadSampler()
    while True: 
        line1 = samfile.readline()
        line2 = samfile.readline()
        if not line2: 
            break
        line1 = SamAlignment(line1)
        line2 = SamAlignment(line2)
        if line1.QNAME != line2.QNAME:
            raise ValueError(
                f"Unpaired read or read with more than one pair "\
                f"encountered. Check your input file. File must "\
                f"be sorted by read name, every read must have "\
                f"a single pair and each pair must have one mapping. "\
                f"{line1.QNAME} {line2.QNAME}."
            )
        if line1.BITS[6]:
            r1 = line1
            r2 = line2
        elif line1.BITS[7]:
            r1 = line2
            r2 = line1
        else:
            raise ValueError(
                f"Could not assign R1 and R2 to reads"\
                f"{line1.QNAME} {line2.QNAME}."
            )
        try:
            # vals is [[start, end, strand]]
            vals = get_paired_blocks(r1,r2)
            ctg_idx = ctg_lut[r1.RNAME]["idx"]
            vals.append(ctg_idx)
            read_sampler.add_read(vals)
        except ValueError as err:
            logging.error("Skipping pair {}".format(err))
        except RuntimeError as err:
            logging.error("Skipping pair {}".format(err))
    return read_sampler

@numba.jit(nopython=True)
def map_read(array, read):
    start, stop, strand = read
    # the below line implements linear scaling with read length
    array[start:stop, strand] += 100.0/(stop-start)

#def sample(read_sampler, n, array):
#    """ Sample reads with replacement from a sampler and map them to an array
#
#    Args:
#    -----
#        read_sampler : ReadSampler object
#            object holding reads to sample
#        n : int
#            number of reads to sample
#        array : 1d np.array
#            
#    Modifies:
#    --------
#        array
#    """
#    sampled_reads = read_sampler.pull_reads(n)
#    fast_sum_coverage(sampled_reads, array)

@numba.jit(nopython=True)
def fast_sum_coverage(reads, array):
    """Map sampled reads to an array using fast jitted function

    Args:
    ------
        reads : list
		    Each element of reads contains [start,end] positions
            for a given aligment
        array : 1d np.array (TO DO: update to 2d, with a
                dimension for +/- strand)
    		array to be populated with coverage

    Modifies:
    -------
        array
    """

    for i in range(reads.shape[0]):
        read = reads[i,:]
        map_read(array, read)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--hdf_file',
        help = "Hdf5 file name which contains chromosome information\
                and to which resulting parser object will be written.",
        required=True,
    )
    parser.add_argument(
        '--global_conf_file',
        help = "Path to global configuration file for this experiment.",
        required = True,
    )
    subparsers = parser.add_subparsers(help='commands', dest='command')

    # parse
    parse_parser = subparsers.add_parser('parse', help="create a sampling\
        object from a sam file")
    
    parse_parser.add_argument('samfile', help="Input samfile, this tool does no\
        filtering and will consider every line in the file. Accepts input\
        from stdin if '-' is specified here.")
    
    parse_parser.add_argument('--paired', action="store_true", help="Consider\
        the sam file as paired. If this flag is specified then the sam file\
        MUST be pre-filtered to have only ONE alignment per pair. Further,\
        there must be NO unpaired reads in the file and the reads must be\
        sorted by read name.")

    # sample
    sample_parser = subparsers.add_parser('sample', help="sample coverage from a\
        sampling object.")
    # removing array_size since we'll get it from the info_pkl instead
    #sample_parser.add_argument('array_size',type=int, help="length of genome")
    sample_parser.add_argument('--num_samples', type=int, default=1,
        help="number of full samples to pull from the sampler, default is 1")
    sample_parser.add_argument('--num_reads', type=int, default=None, 
        help="number of reads to pull for each sample. Default is the size of\
        sampling object.")
    sample_parser.add_argument('--frac_reads', type=float, default=None, 
        help="fraction of reads to pull for each sample. Default is the size of\
        sampling object.")
    sample_parser.add_argument('--identity', action="store_true",
        help="write an array of the actual coverage without sampling, ignores\
        other optional arguments")
    sample_parser.add_argument('--resolution', type=int, default=1,
        help="only keep data for one bp out of this number")

    args = parser.parse_args()
    HDF = args.hdf_file
    conf_dict_global = toml.load(args.global_conf_file)
    res = conf_dict_global["genome"]["resolution"]
    #stranded = hdf_utils.is_stranded(HDF)
    #assert stranded == conf_dict_global["general"]["stranded"], "ERROR: your hdf5 file's 'stranded' attribute and global configuration file's ['general']['stranded'] do not match. Exiting now."

    ctg_lut = hdf_utils.get_ctg_lut(HDF)

    if args.command == "parse":

        if args.samfile == "-":
            print("Reading sam records from stdin")
            f = sys.stdin
        else:
            print(f"Reading sam records from {args.samfile}")
            f = open(args.samfile, mode="r")
        if args.paired:
            try:
                sampler = create_read_list_paired(f, ctg_lut)
            except Exception as e:
                f.close()
                sys.exit(f"ERROR in create_read_list_paired: {e}")
        else:
            try:
                sampler = create_read_list(f, ctg_lut)
            except:
                f.close()
                sys.exit("ERROR in create_read_list")
        f.close()
        sampler.sort_reads()
        hdf_utils.write_dset(HDF, "parser", sampler.reads, np.uint64)

    elif args.command == "sample":
        # Loop over contigs and allocate a separate array for each to hold
        #   coverage calculations. Store them as values in a dictionary with
        #   each contig's index as keys.
        samples_dict = {}
        for ctg_id,ctg_info in ctg_lut.items():
            samples_dict[ctg_info["idx"]] = np.zeros(
                (ctg_info["length"], args.num_samples, 2)
            )

        # Instantialize a ReadSampler obj, and read parser array from hdf5
        sampler = ReadSampler()
        sampler.from_hdf5(HDF)

        if args.identity:
            for ctg_id,ctg_info in ctg_lut.items():
                ctg_reads = sampler.reads[
                    sampler.reads[:,-1] == ctg_info["idx"], 0:3
                ]
                fast_sum_coverage(ctg_reads, samples_dict[ctg_info["idx"]])
                # get every res-th position (rows) of the sampled array,
                # both strands (cols)
                ctg_arr = samples_dict[ctg_info["idx"]][::res,:,:]

                dset_name = "orig"
                hdf_utils.write_dset(
                    HDF,
                    dset_name,
                    ctg_arr,
                    ctg_arr.dtype,
                    group_name = "contigs/{}".format(ctg_id),
                )
            
            bg_outname = HDF.split('.')[0] + "{}_coverage.bedgraph"
            superctg_data = hdf_utils.concatenate_contig_data(
                HDF,
                dset_name,
            )
            hdf_utils.write_bedgraph(
                superctg_data,
                HDF,
                bg_outname.format("both_strand"),
                strand = "both",
            )
            hdf_utils.write_bedgraph(
                superctg_data,
                HDF,
                bg_outname.format("plus_strand"),
                strand = "plus",
            )
            hdf_utils.write_bedgraph(
                superctg_data,
                HDF,
                bg_outname.format("minus_strand"),
                strand = "minus",
            )

        else:

            if args.num_reads:
                num_reads = args.num_reads
            elif args.frac_reads:
                num_reads = int(sampler.total * args.frac_reads)
                print("Total number of reads was {}. Sampling {} reads with replacement.".format(sampler.total, num_reads))
            else:
                num_reads = sampler.total 


            for i in range(args.num_samples):

                sampled_reads = sampler.pull_reads(num_reads)

                progress(
                    i+1,
                    args.num_samples,
                    status='Bootstrap sampling {} times.'.format(args.num_samples)
                )

                for ctg_id,ctg_info in ctg_lut.items():
                    # grab sampled reads of this contig's index
                    ctg_reads = sampled_reads[
                        sampled_reads[:,-1] == ctg_info["idx"], 0:3
                    ]
                    # calculate coverage as each position in contig
                    # Modifies array in dictionary in place
                    fast_sum_coverage(
                        ctg_reads,
                        samples_dict[ctg_info["idx"]][:,i,:]
                    )

            dset_name = "bs"
            for ctg_id,ctg_info in ctg_lut.items():

                ctg_arr = samples_dict[ctg_info["idx"]][::res,:,:]

                hdf_utils.write_dset(
                    HDF,
                    dset_name,
                    ctg_arr,
                    ctg_arr.dtype,
                    group_name = "contigs/{}".format(ctg_id),
                )

            bg_outname = HDF.split('.')[0] + "{}_mean_bootstrapped_coverage.bedgraph"
            superctg_data = hdf_utils.concatenate_contig_data(
                HDF,
                dset_name,
                sample_num = args.num_samples,
            )
            superctg_mean = np.mean(superctg_data, axis=1)
            hdf_utils.write_bedgraph(
                superctg_mean,
                HDF,
                bg_outname.format("both_strands"),
                strand = "both",
            )
            hdf_utils.write_bedgraph(
                superctg_mean,
                HDF,
                bg_outname.format("plus_strands"),
                strand = "plus",
            )
            hdf_utils.write_bedgraph(
                superctg_mean,
                HDF,
                bg_outname.format("minus_strands"),
                strand = "minus",
            )

