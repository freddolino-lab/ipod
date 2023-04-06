#!/usr/bin/python

# script for doing my standardized preprocessing on DNA files
# we do standardized quality trimming and adapter removal

import sys
import subprocess
import tempfile
import os
import toml

# parse the top-level config and get needed general information
conf_dict_global = toml.load(sys.argv[2])

## Set up defined constants that should be universal
proc_opts = conf_dict_global["processing"]

NPROC = 1
if "threads" in proc_opts:
    NPROC = proc_opts["threads"]

PE = conf_dict_global["general"]["paired_reads"]
# set to false by default, and switch only if handle_umi is there, and is set to true.
# this provides backward-compatability
UMI = False
if "handle_umi" in proc_opts:
    UMI = proc_opts["handle_umi"]
if UMI:
    umi_opts = conf_dict_global["umi"]
    UMI_LEN = umi_opts["length"]
    UMI_METHOD = umi_opts["method"]
    UMI_READ = umi_opts["read"]
PROCDIR = proc_opts["processed_direc"]
MAX_ADAPT_N = proc_opts["adapt_max_n"]
TRAILING_JUNK_LEN = proc_opts["trim_trailing_junk_length"]
SLIDE_WIN_LEN = proc_opts["trim_sliding_window_length"]
SLIDE_WIN_QUAL = proc_opts["trim_sliding_window_qual"]
PROCESSED_READ_MINLEN = proc_opts["min_processed_readlen"]
F_READ_SUFFIX = proc_opts["f_paired_read_file_suffix"]
R_READ_SUFFIX = proc_opts["r_paired_read_file_suffix"]
F_UP_READ_SUFFIX = proc_opts["f_unpaired_read_file_suffix"]
R_UP_READ_SUFFIX = proc_opts["r_unpaired_read_file_suffix"]

# base directory for all code for ipod analysis
BINDIR = conf_dict_global["general"]["bindir"]
RAWDIR = conf_dict_global["general"]["rawdir"]
STARTDIR = os.getcwd()

#PREPEND = "{} {{}} {{}} {{}} {{}} {{}}".format(os.path.join(BINDIR, "umi/prepend_umi.sh"))

def concatenate_files(name_wildcard, tmpfile):
    infile_fwd = tmpfile.name
    cmd1 = f"cat {name_wildcard} > {infile_fwd}"
    print("Concatenating files from separate runs")
    print(cmd1)
    res = subprocess.run(cmd1, shell=True, capture_output=True)
    if res.returncode != 0:
        raise Exception(res.stderr)
    
    return infile_fwd

# define some functions that will be used in the rest of the script
def preprocess_file(samp):
    '''Do some initial preprocessing of a gz file,
    including trimming and quality score filtering.
    Don't use this function outside of this script;
    it depends on information defined above.

    Args:
    -----
    samp : dict
        This dict object has all the data listed
        higher up in this file

    Returns:
    --------
    None
        Creates fq.gz files containing processed reads
        Two files for the paired reads (fr, rev)
        Two files for the unpaired reads (fr, rev)
    '''

    infile_1 = samp["ffile"]
    infile_2 = samp["rfile"]
    outprefix = samp["outprefix"]
    PHRED_BASE = samp["phredbase"]
    FWD_ADAP_SEQ = samp["fwd_adapseq"]
    REV_ADAP_SEQ = samp["rev_adapseq"]
    pe = samp["paired_reads"]

    #if infile_1[-3:] == ".gz":
    #    DCPROG = 'zcat'
    #elif infile_1[-4:] == ".bz2":
    #    DCPROG = 'bzcat'
    #else:
    #    raise("Could not determine the decompression program to use")

    #if DCPROG == "bzcat":

    #    in1 = tempfile.NamedTemporaryFile(suffix='.fastq')
    #    infile_fwd = in1.name
    #    cmd1="{} {} > {}".format(DCPROG, infile_1, infile_fwd)
    #    subprocess.call(cmd1,shell=True)

    #    if pe:
    #        in2 = tempfile.NamedTemporaryFile(suffix='.fastq')
    #        infile_rev = in2.name
    #        cmd2="{} {} > {}".format(DCPROG, infile_2, infile_rev)
    #        subprocess.call(cmd2,shell=True)

    #else:

    # get shell output for list of files matching input, since sometimes
    # we have multiple separate starting files which we'd have to concatenate
    files = subprocess.run(f"ls {infile_1}", shell=True, capture_output=True)
    print(files)
    infile_fwd_list = files.stdout.decode().strip().split("\n")
    print(infile_fwd_list)

    # if multiple input fwd files, concatenate to single fastq.qz file
    if len(infile_fwd_list) > 1:
        concat_fwd_infile = tempfile.NamedTemporaryFile(suffix='.fastq.gz')
        infile_fwd = concatenate_files(infile_1, concat_fwd_infile)
    else:
        infile_fwd = infile_fwd_list[0]
    if pe:
        # get shell output for list of files matching input, since sometimes
        # we have multiple separate starting files which we'd have to concatenate
        files = subprocess.run(
            f"ls {infile_2}",
            shell=True,
            capture_output=True,
        )
        infile_rev_list = files.stdout.decode().strip().split("\n")
        if len(infile_rev_list) > 1:
            concat_rev_infile = tempfile.NamedTemporaryFile(
                suffix='.fastq.gz'
            )
            infile_rev = concatenate_files(infile_2, concat_rev_infile)
        else:
            infile_rev = infile_rev_list[0]

    # these get changed below if umi handling is performed
    trim_in_fwd = infile_fwd
    trim_in_rev = infile_rev

    # Deduplicate reads using UMI if that's been selected
    if UMI:

        if UMI_METHOD == "5-prime":

            if UMI_READ == "R2":
                umi_file_1 = infile_rev
                extfile_1 = os.path.join(PROCDIR, outprefix+"_rev_extract.fq.gz")
                trim_in_rev = extfile_1
            elif UMI_READ == "R1":
                umi_file_1 = infile_fwd
                extfile_1 = os.path.join(PROCDIR, outprefix+"_fwd_extract.fq.gz")
                trim_in_fwd = extfile_1

            bc = "N" * UMI_LEN
            if pe:
                if UMI_READ == "R2":
                    umi_file_2 = infile_fwd
                    extfile_2 = os.path.join(PROCDIR, outprefix+"_fwd_extract.fq.gz")
                    trim_in_fwd = extfile_2
                elif UMI_READ == "R1":
                    umi_file_2 = infile_rev
                    extfile_2 = os.path.join(PROCDIR, outprefix+"_rev_extract.fq.gz")
                    trim_in_rev = extfile_2
                ext_cmd = f"umi_tools extract -I {umi_file_1} --bc-pattern={bc} "\
                    f"--read2-in={umi_file_2} --log=umi_extract.log "\
                    f"--read2-out={extfile_2} --stdout={extfile_1}"
            else:
                ext_cmd = f"umi_tools extract --stdin={umi_file_1} --bc-pattern={bc} "\
                    f"--log=umi_extract.log --stdout={extfile_1}"

            res = subprocess.run(ext_cmd, shell=True, capture_output=True)
            if res.returncode != 0:
                print(f"{ext_cmd} returned non-zero exit status:", file=sys.stderr)
                print("stdout:")
                print(res.stdout.decode(), file=sys.stdout)
                print(res.stderr.decode(), file=sys.stderr)
                sys.exit()
            else:
                print(f"{ext_cmd} completed normally with the following output:")
                print(res.stdout.decode(), file=sys.stdout)

##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
## neb method might not need "extract" run if colon can be used. If underscore    
##   must be used then colon would have to be replaced or Peter's demux would need
##   adjusted to put "_" insead of ":"
##################################################################################
##################################################################################
##################################################################################
        elif UMI_METHOD == "NEB":
            # the 11-base UMI is in the read names already, I just need to
            # know whether the ":" can be used, or whether "_" must be used instead
            print()

    cutfile_fwd = os.path.join(PROCDIR, outprefix+"_fwd_cutadap.fq.gz")
    cutfile_rev = os.path.join(PROCDIR, outprefix+"_rev_cutadap.fq.gz")

    # do adapter/quality trimming and write a processed file
    # first clip the adapter sequences from 3-prime ends of reads
    if pe:
        cutadapt_cmd = f"cutadapt -j {NPROC} --quality-base={PHRED_BASE} "\
            f"-a {FWD_ADAP_SEQ} -A {REV_ADAP_SEQ} -n {MAX_ADAPT_N}"\
            f" --match-read-wildcards "\
            f"-o {cutfile_fwd} -p {cutfile_rev} {trim_in_fwd} {trim_in_rev} "
    else:
        cutadapt_cmd = "cutadapt -j {NPROC} --quality-base={PHRED_BASE} "\
            f"-a {FWD_ADAP_SEQ} -n {MAX_ADAPT_N} --match-read-wildcards "\
            f"-o {cutfile_fwd} {trim_in_fwd} "

    cutadapt_cmd += f"> {outprefix}_cutadapt.log 2> {outprefix}_cutadapt.err"
    trim_in_fwd = cutfile_fwd
    trim_in_rev = cutfile_rev

    print("\n{}\n".format(cutadapt_cmd))
    subprocess.call(cutadapt_cmd, shell=True)

    # Next do quality trimming -- trim true crap from the 3' end,
    #   and then look for a sliding window of desired number of bp
    #   with qualities above 15.
    # We avoid doing 5' end trimming so we don't move where the true end
    #  of the read is.
    # Drop the read if we have less than 10 bases after this is done

    trim_fwd_paired = os.path.join( PROCDIR, outprefix + F_READ_SUFFIX )
    trim_fwd_unpaired = os.path.join(
        PROCDIR, outprefix + F_UP_READ_SUFFIX
    )
    trim_rev_paired = os.path.join( PROCDIR, outprefix + R_READ_SUFFIX )
    trim_rev_unpaired = os.path.join(
        PROCDIR, outprefix + R_UP_READ_SUFFIX
    )

    if pe:
        trim_cmd = f"trimmomatic PE -threads {NPROC} -phred{PHRED_BASE} "\
            f"{trim_in_fwd} {trim_in_rev} {trim_fwd_paired} {trim_fwd_unpaired} "\
            f"{trim_rev_paired} {trim_rev_unpaired} "\
            f"TRAILING:{TRAILING_JUNK_LEN} "\
            f"SLIDINGWINDOW:{SLIDE_WIN_LEN}:{SLIDE_WIN_QUAL} "\
            f"MINLEN:{PROCESSED_READ_MINLEN} 2> {outprefix}_trimmomatic.err"
    else:
        trim_cmd = f"trimmomatic SE -threads {NPROC} -phred{PHRED_BASE} "\
            f"{trim_in_fwd} {trim_fwd_paired} "\
            f"TRAILING:{TRAILING_JUNK_LEN} "\
            f"SLIDINGWINDOW:{SLIDE_WIN_LEN}:{SLIDE_WIN_QUAL} "\
            f"MINLEN:{PROCESSED_READ_MINLEN} 2> {outprefix}_trimmomatic.err"
    print(f"\n{trim_cmd}\n")
    subprocess.call(trim_cmd,shell=True)

    #if DCPROG == "bzcat":
    #    in1.close()
    #    in2.close()


conf_dict = toml.load(sys.argv[1]) # this is the condition-level conf file
samp_types = conf_dict["general"]["sample_types"]

for samp_type in samp_types:

    print("Running preprocessing for {} sample.".format(samp_type))

    # gather this sample type's info from conf file
    samp_info = conf_dict[samp_type]
    freads = samp_info["R1_raw_files"]
    rreads = samp_info["R2_raw_files"]
    if "adapter_seqs" in samp_info:
        fwd_adapts = samp_info["adapter_seqs"]
        rev_adapts = samp_info["adapter_seqs"]
    if "R1_adapter_seqs" in samp_info:
        fwd_adapts = samp_info["R1_adapter_seqs"]
    if "R2_adapter_seqs" in samp_info:
        rev_adapts = samp_info["R2_adapter_seqs"]
    rep_names = samp_info["sample_names"]
    sample_direc = samp_info["directory"]

    # move into this sample's directory
    os.chdir(sample_direc)

    # if the rawdir option is not "None" (note a string, not a None object)
    #   then determine whether it's already a symlink. If it's not a symlink
    #   already, create the symlink within the data directory
    if RAWDIR != "None":
        # if it's not a directory, check if it's a symlink
        if not os.path.isdir('raw'):
            if not os.path.islink('raw'):
                os.symlink(conf_dict_global["general"]["rawdir"], "raw")
            # if it is a symlink, remove and re-create it. This ensures path is
            #   still accurate
            else:
                os.remove("raw")
                os.symlink(conf_dict_global["general"]["rawdir"], "raw")

    # set up the needed directories if they are not already present
    if not(os.path.isdir(PROCDIR)):
        os.mkdir(PROCDIR)

    # loop over each replicate's information and do preprocessing
    for i in range(len(rep_names)):
        samp_dict = {
            "ffile": os.path.join("raw", freads[i]),
            "rfile": os.path.join("raw", rreads[i]),
            "fwd_adapseq": fwd_adapts[i],
            "rev_adapseq": rev_adapts[i],
            "phredbase": conf_dict_global["general"]["phredbase"],
            "outprefix": rep_names[i],
            "paired_reads": PE,
        }
        preprocess_file(samp_dict)
    
    # when finished with this sample's preprocessing, move back up
    #  to condition-level directory
    os.chdir(STARTDIR)

