import sys
import os
import h5py
import numpy as np
import anno_tools
import subprocess

# Utilities for working with hdf5 files

def make_ctg_lut_from_bowtie(bowtie_idx):
    '''This function takes a bowtie2 index as its sole argument and
    returns each contig's id and size
    '''

    inspect_cmd = 'bowtie2-inspect -s {} | grep ^Sequence | cut -f 2,3'.format(
        bowtie_idx
    )
    cmd_out = subprocess.check_output(inspect_cmd, shell=True)
    cmd_lines = cmd_out.decode().strip().split('\n')
    cmd_lines.sort()

    ctg_lut = {}
    for ctg_idx,ctg in enumerate(cmd_lines):
        ctg_str,ctg_len = ctg.split('\t')
        # if the contig id field has spaces, split on them and take
        #   first element as ctg_id
        ctg_id = ctg_str.split(' ')[0]
        ctg_len = int(ctg_len)
        ctg_lut[ctg_idx] = {"id":ctg_id, "length":ctg_len}

    return ctg_lut


def set_up_hdf_file(hdf_name, ctg_lut, res, type_lut=None, paired=False):
    '''Write an hdf5 file containing metadata that will be used for organizing
    parser, coverage, and bootstrapped coverage data. This hdf5 file will also
    store the data, which is written to the file by the subprocesses called
    by this script.

    Args:
    -----
    hdf_name : str
        Name of the hdf5 file
    ctg_lut : dict
        Dictionary with desired indices, ids, and length information for
        contigs found in this refernce genome.
    res : int
        Resolution for this experiment.
    '''

    with h5py.File(hdf_name, 'w') as hf:
        # add resolution for this experiment as attribute of the file
        hf.attrs["resolution"] = res
        # create a group called 'contigs'
        ctg_group = hf.create_group('contigs')
        # loop over contigs, gather information for attributes
        for ctg_idx,ctg_info in ctg_lut.items():

            ctg_len = ctg_info["length"]
            ctg_id = ctg_info["id"]
            g_locs = np.arange(0,ctg_len,res)

            this_ctg = ctg_group.create_group(ctg_id)
            this_ctg.attrs["idx"] = ctg_idx
            this_ctg.attrs["length"] = ctg_len
            
            loci = this_ctg.create_dataset(
                "loci",
                g_locs.shape,
                dtype=np.int64,
            )
            loci[...] = g_locs

            if type_lut is not None:
                for samp_type in type_lut.keys():
                    this_samp = this_ctg.create_group(samp_type)
                    if paired:
                        this_samp.create_group('replicates')


def get_resolution(hdf_name):
    with h5py.File(hdf_name, 'r') as hf:
        return hf.attrs["resolution"]


def get_ctg_lut(hdf_name):
    '''Read in hdf5 file and set up contig lookup table.'''

    ctg_lut = {}
    with h5py.File(hdf_name, 'r') as hf:
        # loop over contig groups in contigs
        for ctg_id,ctg_grp in hf['contigs'].items():
            ctg_lut[ctg_id] = {}

            for att_name,att_val in ctg_grp.attrs.items():
                ctg_lut[ctg_id][att_name] = att_val

    return ctg_lut


def load_dset(hdf_name, dset_name, group_name='/'):
    '''Opens hdf file for reading and leads dataset into numpy array.
    '''
    with h5py.File(hdf_name, 'r') as hf:
        dset_path = os.path.join(group_name,dset_name)
        out_arr = hf[dset_path][...]
    return out_arr


def write_dset(hdf_name, dset_name, data_arr, dtype, group_name='/',
               compression='gzip', compression_opts=9, attrs={}):
    '''Opens hdf5 file for appending, creates a dataset, writes data to it,
    and closes the hdf5 file.

    Args:
    -----
    hdf_name : str
        Path to hdf5 file.
    dset_name : str
        Name of dataset that will be created and written to hdf5 file.
    data_arr : np.array
        Data to write to dset.
    dtype : np.dtype
        The type that the data are.
    group_name : str
        Default is '/', or the root of the hierarchy in the hdf5 file.
        If not '/', we assume the group already exists and exit with an
        error message if the desired group doesn't yet exist.
    compression : str
        Type of compression to use.
    compression_opts : int [0-9]
        0 is minimal compression, 9 is maximal compression.
    attrs : dict
        Dictionary of attributes to attach to the written dataset
    '''

    assertion_msg = "compression_opts argument to hdf_utils.write_dset must be an integer from 0 through 9."
    assert ((compression_opts <= 9 and compression_opts >= 0) and isinstance(compression_opts, int)), assertion_msg
    
    with h5py.File(hdf_name, 'a') as hf:
        if not group_name in hf:
            sys.exit("ERROR in hdf_utils.write_dset. Group {} does not exist.".format(group_name))
        grp = hf[group_name]
        # If the dset already exists, we assume it's outdated and replace.
        if dset_name in grp:
            del grp[dset_name]
        dset = grp.create_dataset(
            dset_name,
            data_arr.shape,
            dtype=dtype,
            compression = compression,
            compression_opts = compression_opts,
        )
        dset[...] = data_arr
        for k,v in attrs.items():
            dset.attrs.create(
                name = k,
                data = v,
            )


def calc_supercontig_posnum(hdf_name, spikein_name=None):
    '''Calculates total number of genome positions considered in a 
    supercontig. A supercontig is what we're calling the concatenated
    data from all contigs in a reference sequence.

    If spikein_name is not None (defaults to None), then we skip
    the spike-in contig in determining the number of positions
    in the genome.
    '''

    ctg_lut = get_ctg_lut(hdf_name)

    positions = 0
    with h5py.File(hdf_name, 'r') as hf:

        for ctg_id,ctg_info in ctg_lut.items():
            if ctg_id == spikein_name:
                continue
            # look in contigs/ctg_id group for loci dataset
            loci = hf["contigs/{}/loci".format(ctg_id)].shape[0]
            positions += loci

    return positions

def calc_spikein_posnum(hdf_name, spikein_name):
    
    ctg_lut = get_ctg_lut(hdf_name)

    positions = 0
    with h5py.File(hdf_name, 'r') as hf:

        for ctg_id,ctg_info in ctg_lut.items():
            # only count up positions is we're looking at the spike-in
            if ctg_id == spikein_name:
                loci = hf["contigs/{}/loci".format(ctg_id)].shape[0]
                positions += loci

    return positions


def concatenate_contig_data(hdf_name, dset_basename="orig", sample_num=1,
                            positions=None, spikein_name=None):
    '''Loops over contigs in ctg_lut, reads orig data from hdf5 file for each,
    and appends values to long super-contig. Returns the super-contig's vals.

    If spikein_name is not None, then this indicates that one contig in the
    reference sequence was for spike-in data to align. In that case, peel
    off the spike-in counts to a separate array so we can do something else
    with that information.
    '''

    ctg_lut = get_ctg_lut(hdf_name)

    if positions is None:
        positions = calc_supercontig_posnum(hdf_name, spikein_name)
    sample_arr = np.zeros((positions,sample_num))

    prior_stop = 0

    with h5py.File(hdf_name, 'r') as hf:
        for ctg_id,ctg_info in ctg_lut.items():
            dset_name = "contigs/{}/{}".format(ctg_id, dset_basename)
            these_vals = hf[dset_name]
            # only grab data if the contig is not the spike-in sequence
            if ctg_id != spikein_name:
                # look in contigs/ctg_id group for loci dataset
                ctg_positions = hf["contigs/{}/loci".format(ctg_id)].shape[0]
                current_stop = prior_stop + ctg_positions
                sample_arr[prior_stop:current_stop,:] = these_vals
                prior_stop = current_stop

    return sample_arr

def get_spikein_data(hdf_name, spikein_name, dset_basename="orig"):
    '''Gets spike-in data from the hdf5 file.
    
    Args:
    -----
    hdf_name : str
        Name of the hdf5 file containing data.
    spikein_name : str
        Name of the "chromosome" in your reference file representing
        the spike-in sequence.
    dest_basename : str
        Dataset name to grab data from. Default is "orig".
    sample_num : int
        The number of samples drawn at bootstrapping step.
    '''

    ctg_lut = get_ctg_lut(hdf_name)

    dset_name = "contigs/{}/{}".format(spikein_name, dset_basename)
    with h5py.File(hdf_name, 'r') as hf:
        spikein_arr = hf[dset_name][...]

    return spikein_arr

def decatenate_supercontig_data(hdf_name, superctg_arr, ctg_lut):
    '''Splits data in superctg_arr into each contig's values.

    Args:
    -----
    hdf_name : str
        Path to hdf file.
    superctg_arr : np.array
        Array containing data for all contigs appended to each other
    ctg_lut : dict
        Dictionary mapping contig names to other information about each contig.
    
    Returns:
    --------
    ctg_data : dict
        Dictionary to map contig names to the appropriate data.
    '''
    ctg_data = {}
    start_idx = 0
    for ctg_id,ctg_info in ctg_lut.items():
        with h5py.File(hdf_name, 'r') as hf:
            ctg_positions = hf["contigs/{}/loci".format(ctg_id)].shape[0]

        end_idx = ctg_positions + start_idx
        if superctg_arr.ndim == 1:
            superctg_arr = np.expand_dims(superctg_arr, -1)
        ctg_data[ctg_id] = superctg_arr[start_idx:end_idx,:]
        start_idx = end_idx

    return ctg_data
    

def decatenate_and_write_supercontig_data(hdf_name, superctg_arr,
    dset_name, dtype=np.float64, grp_fmt_str="contigs/{}", attrs={}):
    '''Splits data in a supercontig array into each original contig's values.
    '''
    ctg_lut = get_ctg_lut(hdf_name)
    ctg_data = decatenate_supercontig_data(hdf_name, superctg_arr, ctg_lut)

    for ctg_id,ctg_vals in ctg_data.items():
        grp = grp_fmt_str.format(ctg_id)
        write_dset(
            hdf_name,
            dset_name,
            ctg_vals,
            dtype,
            group_name = grp,
            attrs = attrs,
        )


def create_wig_record(superctg_arr, hdf_name):
    '''Decatenates data in superctg_arr into its appropriate contigs
    and returns a WigData object.

    Args:
    -----
    superctg_arr : np.array
        Array with data for all contigs appended to each other.
    hdf_name : str
        Path to hdf file.

    Returns:
    --------
    wig_record : anno_tools.WigData
    '''

    ctg_lut = get_ctg_lut(hdf_name)
    ctg_data = decatenate_supercontig_data(hdf_name, superctg_arr, ctg_lut)
    resolution = get_resolution(hdf_name)

    wig_record = anno_tools.WigData()
    #################################################################
    # NOTE: add data_origin as cfg to ensure it's always up-to-date!!
    #################################################################
    for ctg_id,ctg_vals in ctg_data.items():
        with h5py.File(hdf_name, 'r') as hf:
            ctg_locs = hf["contigs/{}/loci".format(ctg_id)]

            for i,loc in enumerate(ctg_locs):
                start = loc
                end = loc + resolution 
                wig_record.addline(
                    chrom_name = ctg_id,
                    start = start,
                    end = end,
                    score = ctg_vals[i][0],
                )

    return wig_record


def write_bedgraph(superctg_arr, hdf_name, out_fname, spikein_name=None):
    '''Decatenates data in superctg_arr into its appropriate contigs
    and returns a FastBEDGraphData object.

    Args:
    -----
    superctg_arr : np.array
        Array with data for all contigs appended to each other.
    hdf_name : str
        Path to hdf file.
    out_fname: str
        Path to output bedgraph file.
    spikein_name : str
        Name of spike-in chromosome.

    Returns:
    --------
    bed_record : anno_tools.FastBEDGraphData
    '''

    ctg_lut = get_ctg_lut(hdf_name)
    ctg_data = decatenate_supercontig_data(
        hdf_name,
        superctg_arr,
        ctg_lut,
    )
    resolution = get_resolution(hdf_name)

    with open(out_fname, 'w') as outfile:

        for ctg_id,ctg_vals in ctg_data.items():
            if ctg_id == spikein_name:
                continue
            with h5py.File(hdf_name, 'r') as hf:
                ctg_locs = hf["contigs/{}/loci".format(ctg_id)][...]

            ends = ctg_locs + resolution
            scores = ctg_vals[:,0]

            for i in range(len(ctg_locs)):
                outfile.write(
                    "{}\t{}\t{}\t{}\n".format(
                        ctg_id,
                        ctg_locs[i],
                        ends[i],
                        scores[i],
                    )
                )


def create_bedgraph_record(superctg_arr, hdf_name):
    '''Decatenates data in superctg_arr into its appropriate contigs
    and returns a BEDGraphData object.

    Args:
    -----
    superctg_arr : np.array
        Array with data for all contigs appended to each other.
    hdf_name : str
        Path to hdf file.

    Returns:
    --------
    bed_record : anno_tools.BEDGraphData
    '''

    ctg_lut = get_ctg_lut(hdf_name)
    ctg_data = decatenate_supercontig_data(hdf_name, superctg_arr, ctg_lut)
    resolution = get_resolution(hdf_name)

    bed_record = anno_tools.BEDGraphData()
    for ctg_id,ctg_vals in ctg_data.items():
        with h5py.File(hdf_name, 'r') as hf:
            ctg_locs = hf["contigs/{}/loci".format(ctg_id)]

            for i,loc in enumerate(ctg_locs):
                start = loc
                end = loc + resolution 
                bed_record.addline(
                    chrom_name = ctg_id,
                    start = start,
                    end = end,
                    score = ctg_vals[i][0],
                )

    return bed_record


def create_group(hdf_name, group_name):
    with h5py.File(hdf_name, 'a') as hf:
        if not group_name in hf:
            hf.create_group(group_name)
