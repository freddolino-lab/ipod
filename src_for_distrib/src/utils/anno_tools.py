#!usr/bin/python

import numpy as np
from scipy.spatial import distance
import operator
import tempfile
import subprocess
import os
import sys

# additional functions for splitting the last field of a gff entry
def newSplit(value):
        import shlex
        lex = shlex.shlex(value)
        lex.quotes = '"'
        lex.whitespace_split = True
        lex.commenters = ''
        return list(lex)

def make_comment_dict(gff_entry):
        gff_entry.comment_dict = {}
        keyvalue = gff_entry.comments.split(";")
        for pair in keyvalue:
                key_value= newSplit(pair)
                key = key_value[0]
                if len(key_value) < 2:
                        value = ""
                else:
                        value = "".join(key_value[1:])
                gff_entry.comment_dict[key] = value

class NarrowPeakEntry:
    '''A container class for the equivalent of a narrowpeak file line.

    Attributes:
    -----------
    chrom_name : str
        Name of the chromosome on which peak is found 
    start : int
        Zero-indexed start position of peak region
    end : int
        Zero-indexed end position of peak region
    name : str
        Name of peak, usually not applicable, in which case '.' is used.
    display : int (0-1000)
        Denotes how dark the peak should display in a genome browser.
        Called 'score' on UCSC explanation of NarrowPeak format, but we use
        the term 'score' to mean the median enrichment score within the peak.
    strand : str
        '+' if on plus strand, '-' if on minus strand, '.' if not stranded.
    score : float
        Median enrichment score within the locus specified by start-end.
    pval : float
        p-value for peak, set to -1 if not applicable.
    qval : float
        fdr-corrected significance value. Set to -1 if not applicable.
    peak : float
        The single position where the peak could be considered to exist.
        Usually defined as either the position where the max peak height is
        realized, or could simply be the center position of the peak.
    '''

    VANILLA_FORMAT_STRING = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
    IDR_FORMAT_STRING = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
   
    def __init__(self, line=None):

        if line is not None:
            self.parse_narrowpeak_line(line)
        else:
            self.chrom_name = '' 
            self.start = 0 # zero-indexed
            self.end = 0
            self.name = '.'
            self.display = 0
            self.strand = "."
            self.score = 0
            self.pval = -1
            self.qval = -1
            self.peak = -1 # zero-indexed position of peak
            self.local_idr = -1
            self.global_idr = -1
            self.repa_start = -1
            self.repa_end = -1
            self.repa_score = -1
            self.repa_peak = -1
            self.repb_start = -1
            self.repb_end = -1
            self.repb_score = -1
            self.repb_peak = -1

    def copy(self):
        new = NarrowPeakEntry()
        attr_dict = self.__dict__
        for k,v in attr_dict.items():
            setattr(new, k, v)
        return new

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __len__(self):
        return len(self.__dict__)

    def parse_narrowpeak_line(self, line):
        """
        Set this entry's values to those of a line from a gff file
        """
        e = line.rstrip().split("\t")
        e[1] = int(e[1])
        e[2] = int(e[2])
        e[3] = str(e[3])
        e[4] = float(e[4])
        e[6] = float(e[6])
        e[7] = float(e[7])
        e[8] = float(e[8])
        e[9] = float(e[9])
        if len(e) == 10:
            (
                self.chrom_name,
                self.start,
                self.end,
                self.name,
                self.display,
                self.strand,
                self.score,
                self.pval,
                self.qval,
                self.peak
            ) = e
            self.repb_peak = None
            
        elif len(e) == 20:
            e[10] = float(e[10])
            e[11] = float(e[11])
            e[12] = int(e[12])
            e[13] = int(e[13])
            e[14] = float(e[14])
            e[15] = float(e[15])
            e[16] = int(e[16])
            e[17] = int(e[17])
            e[18] = float(e[18])
            e[19] = float(e[19])
            (
                self.chrom_name,
                self.start,
                self.end,
                self.name,
                self.display,
                self.strand,
                self.score,
                self.pval,
                self.qval,
                self.peak,
                self.local_idr,
                self.global_idr,
                self.repa_start,
                self.repa_end,
                self.repa_score,
                self.repa_peak,
                self.repb_start,
                self.repb_end,
                self.repb_score,
                self.repb_peak,
            ) = e

    def filter(self, attr, val):
        return getattr(self, attr) == val

    def __repr__(self):
        """
        Return a formatted bedgraph line, which can be used to reconstitute the object or be written directly to a bedgraph file
        """
        if self.repb_peak is None:
            return NarrowPeakEntry.VANILLA_FORMAT_STRING.format(
                self.chrom_name,
                self.start,
                self.end,
                self.name,
                self.display,
                self.strand,
                self.score,
                self.pval,
                self.qval,
                self.peak
            )
        else:
            return NarrowPeakEntry.IDR_FORMAT_STRING.format(
                self.chrom_name,
                self.start,
                self.end,
                self.name,
                self.display,
                self.strand,
                self.score,
                self.pval,
                self.qval,
                self.peak,
                self.local_idr,
                self.global_idr,
                self.repa_start,
                self.repa_end,
                self.repa_score,
                self.repa_peak,
                self.repb_start,
                self.repb_end,
                self.repb_score,
                self.repb_peak
            )


class BEDEntry:
    '''A simple container class for the equivalent of a bed6 file line.
    '''
    # contig_name start end sitename score
    FORMAT_STRING = "{}\t{}\t{}\t{}\t{}\t{}"

    def __init__(self, line=None):

        if line is not None:
            # set defaults, then reassign depending on line
            self.chrom_name = ""
            self.start = 0
            self.end = 0
            self.name = ""
            self.score = 0
            self.strand = "."
            self.parse_bed_line(line)
        else:
            self.chrom_name = ""
            self.start = 0
            self.end = 0
            self.name = ""
            self.score = 0
            self.strand = "."

    def parse_bed_line(self, line):
        """
        Set this entry's values to those of a line from a bed6 file
        """

        datarray = line.rstrip().split("\t")
        type_mapper = {
            'chrom_name': str,
            'start': int,
            'end': int,
            'name': str,
            'score': float,
            'strand': str,
        }
        attr_list = ['chrom_name', 'start', 'end', 'name', 'score', 'strand']
        while len(datarray) > 0:
            field_info = datarray.pop(0)
            field = attr_list.pop(0)
            field_type = type_mapper[field]
            setattr(self, field, field_type(field_info))

    def filter(self, attr, val):
        return getattr(self, attr) == val

    def __repr__(self):
        """
        Return a formatted bed line, which can be used to reconstitute the object or be written directly to a bed file
        """
        return BEDEntry.FORMAT_STRING.format(
            self.chrom_name,
            self.start,
            self.end,
            self.name,
            self.score,
            self.strand,
        )


class BEDGraphEntry:
    '''A simple container class for the equivalent of a bedgraph file line.
    '''
    # contig_name start end sitename score
    FORMAT_STRING = "{}\t{}\t{}\t{}"

    def __init__(self, line=None):

        if line is not None:
            self.parse_bedgraph_line(line)
        else:
            self.chrom_name = ""
            self.start = 0
            self.end = 0
            self.score = 0

    def parse_bedgraph_line(self, line):
        """
        Set this entry's values to those of a line from a bedgraph file
        """

        datarray = line.rstrip().split("\t")
        self.chrom_name = datarray[0]
        self.start = int(datarray[1])
        self.end = int(datarray[2])
        self.score = float(datarray[3])

    def filter(self, attr, val):
        return getattr(self, attr) == val

    def __repr__(self):
        """
        Return a formatted bedgraph line, which can be used to reconstitute the object or be written directly to a bedgraph file
        """
        return BEDGraphEntry.FORMAT_STRING.format(
            self.chrom_name,
            self.start,
            self.end,
            self.score
        )


class WigEntry:
    '''A simple container class for the equivalent of a wiggle file line.
    '''
    # contig_name start end sitename score
    FORMAT_STRING = "{}\t{}"

    def __init__(self, line=None):

        if line is not None:
            self.parse_wig_line(line)
        else:
            self.chrom_name = ""
            self.start = 0
            self.end = 0
            self.length = self.end - self.start
            self.score = 0

    #def parse_wig_line(self, line):
    #    """
    #    Set this entry's values to those of a line from a gff file
    #    """
    #
    #    datarray = line.rstrip().split("\t")
    #    self.chrom_name = datarray[0]
    #    self.start = datarray[1]
    #    self.end = datarray[2]
    #    self.score = datarray[4]

    def __repr__(self):
        """
        Return a formatted wig line, which can be used to reconstitute the object or be written directly to a bedgraph file
        """
        return WigEntry.FORMAT_STRING.format(
            self.start,
            self.score
        )


class AnnotationData:
    '''Class with some utilities for manipulating genome annotation files
    '''

    def __init__(self):
        self.data = []
        self.index = 0
        self.fname = None

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == len(self.data):
            self.index = 0
            raise StopIteration
        self.index = self.index + 1
        return self.data[self.index-1]

    def __getitem__(self, index):
        return self.data[index]

    def __len__(self):
        return len(self.data)

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def clear_db(self):
        self.data = []

    def sort(self):
        """
        sort entries based on contig and starting position
        """
        self.data = sorted(self.data, key=lambda x: (x.chrom_name, x.start))
  
    def find_entry(self, findfunc, findall=False):
        """
        Return the line or lines for which findfunc is true when given an
        annotation item.

        If findall is false, only the first such entry is returned

        Findfunc should take a gffline object as its only argument
        """
        matches = filter(findfunc, self.data)

        if len(matches) == 0:
            return []

        if (findall):
            return matches
        else:
            return matches[0]

    def add_entry(self, new_entry):
        # Add an externally constructed Entry object
        self.data.append(new_entry)

    def fetch_array(self, attr='score'):
        vals = [ getattr(entry, attr) for entry in self.data ]
        return np.asarray(vals)

    def ctg_names(self):
        return set([ entry.chrom_name for entry in self])

    def fetch_complement_bed_data(self, contig_lut, filter_chrs=["None"]):
        """Return the complement of the ranges in self

        Args:
        -----
        contig_lut: dict
            Lookup table with contig names as keys and values as a
            dictionary of information about each contig. One of the keys
            to this dictionary of information must be "length", the value
            of which is an integer indicating the length of the contig.
        filter_chrs: list
            Names of chromosomes to omit from the .genome file
        """

        bed_tmp = tempfile.TemporaryDirectory()
        bed_fname = os.path.join(bed_tmp.name, "bedfile.bed")
        genome_fname = os.path.join(bed_tmp.name, "genome.genome")
    
        with open(genome_fname, 'w') as genomef:
            for ctg_idx,ctg_info in contig_lut.items():
                # skip the filtered chromosomes
                if ctg_info['id'] in filter_chrs:
                    continue
                genomef.write("{}\t{}\n".format(ctg_info['id'], ctg_info["length"]))

        #subprocess.run("cat {}".format(genome_fname), shell=True)

        complement_cmd = "bedtools complement -i {} -g {} > {}".format(
                self.fname,
                genome_fname,
                bed_fname,
            )
        retcode = subprocess.call(complement_cmd, shell=True)

        if retcode != 0:
            print("ERROR in writing complement bed file")
            print("Command attempted was:\n{}".format(complement_cmd))
            sys.exit()

        bed = BEDData()
        bed.parse_bed_file(bed_fname)
        bed.sort()

        # remove temp dir and all contents
        bed_tmp.cleanup()

        return(bed)
                
    def write_file(self):
        """
        Write the [sorted] contents of the data to a file
        """
        self.sort()
        with open(self.fname, "w") as ostr:
            for line in self:
                ostr.write("{}\n".format(line))


class NarrowPeakData(AnnotationData):
    """
    Class for storing and manipulating narrowpeak data
    """

    # super().__init__() keeps all parent attributes and methods
    def __init__(self):
        super().__init__()

    def parse_narrowpeak_file(self, filename, clear=True):
        """
        Parse a narrowpeak file and store the lines in self.data

        The current contents of this object are overwritten iff clear 
        """

        self.fname = filename
        if clear:
            self.clear_db()

        with open(self.fname, "r") as instr:
            for line in instr:
                if line.startswith("#"):
                    continue

                newline = NarrowPeakEntry(line)
                self.data.append(newline)

    
    def addline(self, chrom_name, start, end, score, name='.', display=0, strand='.', pval=-1, qval=-1, peak=-1, local_idr=None, global_idr=None, repa_start=None, repa_end=None, repa_score=None, repa_peak=None, repb_start=None, repb_end=None, repb_score=None, repb_peak=None):
        """
        Add a line with the given data
        """

        newobj = NarrowPeakEntry()
        newobj.chrom_name = chrom_name
        newobj.start = int(start)
        newobj.end = int(end)
        newobj.display = int(display)
        newobj.name = name
        newobj.strand = strand
        newobj.score = score
        newobj.pval = pval
        newobj.qval = qval
        newobj.peak = peak
        newobj.local_idr = local_idr
        newobj.global_idr = global_idr
        newobj.repa_start = repa_start
        newobj.repa_end = repa_end
        newobj.repa_score = repa_score
        newobj.repa_peak = repa_peak
        newobj.repb_start = repb_start
        newobj.repb_end = repb_end
        newobj.repb_score = repb_score
        newobj.repb_peak = repb_peak

        self.data.append(newobj)

    def filter(self, attr, relate, val):
        '''Filters records in self where val relates to the value in attr
        by the relationship defined by relate.

        Args:
        -----
        attr : str
            Filed name to filter data by
        relate : operator
            An operator function to define how to relate the filter_val to
            the value found in filter_field. Default is operator.ge, so
            by default we count regions where the value in filter_field
            is greather than or equal to filter_val. Valid functions
            include: operator.lt, operator.le, operator.eq, operator.ne,
            operator.ge, operator.gt.
        val : float or str
            Value to compare to the data in attr.
        '''

        ctg_np = NarrowPeakData()
        for record in self:
            if relate(getattr(record, attr), val):
                ctg_np.add_entry(record)
        return ctg_np

    def get_bool_arr_dict(self, ctg_len_dict):

        bool_dict = {}
        for ctg_id,ctg_len in ctg_len_dict.items():
            bool_dict[ctg_id] = np.zeros(ctg_len, dtype='bool')
            ctg_np = self.filter("chrom_name", operator.eq, ctg_id)
            for record in ctg_np:
                for i in range(record.start, record.end):
                    bool_dict[ctg_id][i] = True

        return bool_dict

    def calc_jaccard(self, other, ctg_len_dict):

        this_bool_dict = self.get_bool_arr_dict(ctg_len_dict)
        that_bool_dict = other.get_bool_arr_dict(ctg_len_dict)

        for i,(ctg_id,ctg_len) in enumerate(ctg_len_dict.items()):
            
            if i == 0:
                this_bool = this_bool_dict[ctg_id]
                that_bool = that_bool_dict[ctg_id]
            else:
                this_bool = np.append(this_bool, this_bool_dict[ctg_id])
                that_bool = np.append(that_bool, that_bool_dict[ctg_id])

        return distance.jaccard(this_bool, that_bool)


class WigData(AnnotationData):
    """
    Class for storing and manipulating bedgraph data
    """

    # super().__init__() keeps all parent attributes and methods
    def __init__(self):
        super().__init__()

    def addline(self, chrom_name, start, end, score):
        """
        Add a line with the given data
        """

        newobj = WigEntry()
        newobj.chrom_name = chrom_name
        newobj.start = start
        newobj.end = end
        newobj.score = score
        newobj.length = end - start

        self.data.append(newobj)

    def write_file(self, filename):
        self.sort()
        lengths = [rec.length for rec in self.data]
        length_count = len(set(lengths))
        if length_count == 1:
            fixed_step = True
        else:
            fixed_step = False

        with open(filename, "w") as ostr:
            ostr.write("track type=wiggle_0\n")
            prior_chrom = ''
            for line in self:
                if line.chrom_name != prior_chrom:
                    if fixed_step:
                        ostr.write(
                            "fixedStep chrom={} start=0 step={}\n".format(
                                line.chrom_name,lengths[0]
                            )
                        )
                    else:
                        ################## implement span in future
                        ostr.write("variableStep chrom={}\n".format(line.chrom_name))
                if fixed_step:
                    ostr.write("{}\n".format(line.score))
                else:
                    ostr.write("{}\t{}\n".format(line.start,line.score))
                prior_chrom = line.chrom_name


class BEDData(AnnotationData):
    """
    Class for storing and manipulating bedgraph data
    """

    # super().__init__() keeps all parent attributes and methods
    def __init__(self):
        super().__init__()

    def parse_bed_file(self, filename, clear=True):
        """
        Parse a bed file and store the lines in self.data

        The current contents of this object are overwritten iff clear 
        """

        self.fname = filename
        if clear:
            self.clear_db()

        with open(self.fname, "r") as instr:
            for line in instr:
                if line.startswith("#"):
                    continue

                newline = BEDEntry(line)
                self.data.append(newline)
        self.sort()

    def addline(self, chrom_name, start, end, name, score, strand):
        """
        Add a line with the given data
        """

        newobj = BEDGraphEntry()
        newobj.chrom_name = chrom_name
        newobj.start = int(start)
        newobj.end = int(end)
        newobj.name = name
        newobj.score = float(score)
        newobj.strand = strand

        self.data.append(newobj)


class BEDGraphData(AnnotationData):
    """
    Class for storing and manipulating bedgraph data
    """

    # super().__init__() keeps all parent attributes and methods
    def __init__(self):
        super().__init__()

    def filter(self, attr, relate, val):
        '''Filters records in self where val relates to the value in attr
        by the relationship defined by relate.

        Args:
        -----
        attr : str
            Filed name to filter data by
        relate : operator
            An operator function to define how to relate the filter_val to
            the value found in filter_field. Default is operator.ge, so
            by default we count regions where the value in filter_field
            is greather than or equal to filter_val. Valid functions
            include: operator.lt, operator.le, operator.eq, operator.ne,
            operator.ge, operator.gt.
        val : float or str
            Value to compare to the data in attr.
        '''

        bd = BEDGraphData()
        for record in self:
            if relate(getattr(record, attr), val):
                bd.add_entry(record)
        return bd

    def parse_bedgraph_file(self, filename, clear=True):
        """
        Parse a bed file and store the lines in self.data

        The current contents of this object are overwritten iff clear 
        """

        self.fname = filename
        if clear:
            self.clear_db()

        with open(self.fname, "r") as instr:
            for line in instr:
                if line.startswith("#"):
                    continue

                newline = BEDGraphEntry(line)
                self.data.append(newline)
        self.sort()

    def addline(self, chrom_name, start, end, score):
        """
        Add a line with the given data
        """

        newobj = BEDGraphEntry()
        newobj.chrom_name = chrom_name
        newobj.start = int(start)
        newobj.end = int(end)
        newobj.score = float(score)

        self.data.append(newobj)


class GffEntry:
    """
    A simple container class for the equivalent of a gff file line
    """

    FORMAT_STRING = "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}"

    def __init__(self, line=None):

        if line is not None:
            self.parse_gff_line(line)
        else:
            self.chrom_name = ""
            self.data_origin = ""
            self.site_type = ""
            self.start = 0
            self.end = 0
            self.direction = "+"
            self.comments = ""


    def parse_gff_line(self, line):
        """
        Set this entry's values to those of a line from a gff file
        """

        datarray = line.rstrip().split("\t")
        self.chrom_name = datarray[0]
        self.data_origin = datarray[1]
        self.site_type = datarray[2]
        self.start = int(datarray[3])
        self.end = int(datarray[4])
        self.direction = datarray[6]
        self.comments = " ".join(datarray[8:])
        self.comments = self.comments.replace("\t", " ")
        self.comments = self.comments.replace("\n", "")

    def __repr__(self):
        """
        Return a formatted gff line, which can be used to reconstitute the object or be written directly to a gff file
        """
        return GffEntry.FORMAT_STRING.format(self.chrom_name, self.data_origin, self.site_type, self.start, self.end, self.direction, self.comments)


class GffData:
    """
    Class for storing and manipulating gff data
    """
    # super().__init__() keeps all parent attributes and methods
    def __init__(self):
        super().__init__()

    def parse_gff_file(self,filename, clear=True):
        """
        Parse a gff file and store the lines in self.data

        The current contents of this object are overwritten iff clear 
        """

        if (clear):
            self.clear_db()

        with open(filename, "r") as instr:
            for line in instr:
                if line[0] == "#":
                        continue

                newline = GffEntry(line)
                self.data.append(newline)

    def addline(self, chrom_name, data_origin, site_type,
                start, end, direction, comments):
        """
        Add a line with the given data
        """

        newobj = GffEntry()
        newobj.chrom_name = chrom_name
        newobj.data_origin = data_origin
        newobj.site_type = site_type
        newobj.start = start
        newobj.end = end
        #newobj.score = score
        newobj.direction = direction
        newobj.comments = comments

        self.data.append(newobj)


