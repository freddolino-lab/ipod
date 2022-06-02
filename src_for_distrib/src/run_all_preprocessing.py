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
PE = conf_dict_global["general"]["paired_reads"]
UMI = proc_opts["handle_umi"]
if UMI:
    umi_opts = conf_dict_global["umi"]
    UMI_LEN = umi_opts["length"]
    PARDRE_L = umi_opts["pardre_l"]
    PARDRE_C = umi_opts["pardre_c"]
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

aln_opts = conf_dict_global["alignment"]
NPROC = aln_opts["align_threads"]
MIN_FRAG_LEN = aln_opts["min_fragment_length"]
MAX_FRAG_LEN = aln_opts["max_fragment_length"]

# base directory for all code for ipod analysis
BINDIR = conf_dict_global["general"]["bindir"]
RAWDIR = conf_dict_global["general"]["rawdir"]
STARTDIR = os.getcwd()

PARDIR = os.environ["PARDIR"]
PARBIN = os.path.join(PARDIR, "ParDRe")

if not os.path.isfile(PARBIN):
    raise Exception(
        f"ERROR: ParDRe binary, {PARBIN}, does not exist. "\
        f"Make sure ParDRe is compiled and the PARDIR environment "\
        f"variable is set to the location containing the ParDRe binary."
    )
    sys.exit()

PARDRE = "{} -i {{}} -p {{}} -z \
    -o {{}} -r {{}} \
    -l {{}} -c {{}} \
    > {{}}_pardre.log 2> {{}}_pardre.err".format(PARBIN)
PREPEND = "{} {{}} {{}} {{}} {{}}".format(os.path.join(BINDIR, "umi/prepend_umi.sh"))

def concatenate_files(name_list, tmpfile):
    infile_fwd = tmpfile.name
    cmd1 = f"cat {' '.join(name_list)} > {infile_fwd}"
    print("Concatenating files from separate runs")
    print(cmd1)
    res = subprocess.run(cmd1, shell=True, check=True, capture_output=True)
    
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
    ADAP_SEQ = samp["adapseq"]
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
    infile_fwd = files.stdout.decode().strip().split("\n")

    # if multiple input fwd files, concatenate to single fastq.qz file
    if len(infile_fwd) > 1:
        concat_fwd_infile = tempfile.NamedTemporaryFile(suffix='.fastq.gz')
        infile_fwd = concatenate_files(infile_fwd, concat_fwd_infile)
    else:
        infile_fwd = infile_fwd[0]
    if pe:
        # get shell output for list of files matching input, since sometimes
        # we have multiple separate starting files which we'd have to concatenate
        files = subprocess.run(f"ls {infile_2}", shell=True, capture_output=True)
        infile_rev = files.stdout.decode().strip().split("\n")
        if len(infile_rev) > 1:
            concat_rev_infile = tempfile.NamedTemporaryFile(suffix='.fastq.gz')
            infile_rev = concatenate_files(infile_rev, concat_rev_infile)
        else:
            infile_rev = infile_rev[0]

    cutfile_fwd = os.path.join(PROCDIR, outprefix+"_fwd_cutadap.fq.gz")
    cutfile_rev = os.path.join(PROCDIR, outprefix+"_rev_cutadap.fq.gz")

    # do some quality trimming and write a processed file
    # first clip the adapter sequences from 3-prime ends of reads
    if pe:
        cutadapt_cmd = f"cutadapt --quality-base={PHRED_BASE} "\
            f"-a {ADAP_SEQ} -A {ADAP_SEQ} -n {MAX_ADAPT_N} --match-read-wildcards "\
            f"-o {cutfile_fwd} -p {cutfile_rev} {infile_fwd} {infile_rev} "
    else:
        cutadapt_cmd = "cutadapt --quality-base={PHRED_BASE} "\
            f"-a {ADAP_SEQ} -n {MAX_ADAPT_N} --match-read-wildcards "\
            f"-o {cutfile_fwd} {infile_fwd} "
    # if we're handling UMI's we need to filter by minimum read length
    if UMI:
        cutadapt_cmd += f"-m {PARDRE_L} "

    cutadapt_cmd += f"> {outprefix}_cutadapt.log 2> {outprefix}_cutadapt.err"
    # these get changed below if umi handling is performed
    trim_in_fwd = cutfile_fwd
    trim_in_rev = cutfile_rev

    print("\n{}\n".format(cutadapt_cmd))
    subprocess.call(cutadapt_cmd, shell=True)

    # Next, deduplicate reads using UMI if that's been selected
    if UMI:
        # instantiate pardre infile names,
        # which will be changed later if UMI_METHOD is "NEB"
        pardre_fwd_infile = cutfile_fwd
        pardre_rev_infile = cutfile_rev

        # pardre output file names
        dedupfile_fwd = os.path.join(PROCDIR, outprefix+"_fwd_dedup.fq.gz")
        dedupfile_rev = os.path.join(PROCDIR, outprefix+"_rev_dedup.fq.gz")

        # deduped, de-umi'ed output file names
        trim_in_fwd = os.path.join(PROCDIR, outprefix+"_fwd_cut_dedup.fq.gz")
        trim_in_rev = os.path.join(PROCDIR, outprefix+"_rev_cut_dedup.fq.gz")

        if UMI_READ == "R2":
            cut_U = UMI_LEN
            cut_u = 0
        elif UMI_READ == "R1":
            cut_U = 0
            cut_u = UMI_LEN

        if UMI_METHOD == "NEB":
            # to prepare reads for ParDRe,
            # get the 11-base UMI from the 3-prime end of I1 read,
            # place them on the 5-prime end of either the R2 or the R1 read
            index_fname = infile_fwd.replace("_R1_", "_I1_")
            prep_out_name = os.path.join(PROCDIR, "umi_read.fq.gz")
            src_dir = os.path.join(BINDIR, "umi")

            # prep_cmd below will output umi_read.fastq.gz
            if UMI_READ == "R2":
                read_fname = cutfile_rev
                # clobber pardre_rev_infile now
                pardre_rev_infile = prep_out_name
            elif UMI_READ == "R1":
                read_fname = cutfile_fwd
                # clobber pardre_fwd_infile now
                pardre_fwd_infile = prep_out_name

            prep_cmd = PREPEND.format(index_fname, read_fname, src_dir, PROCDIR)
            print(
                f"Prepending UMI from I1 reads to 5-prime "\
                f"end of {UMI_READ} reads:\n"\
                f"{prep_cmd}"
            )

            res = subprocess.run(prep_cmd, shell=True, capture_output=True)
            if res.returncode != 0:
                print(f"{prep_cmd} returned non-zero exit status:", file=sys.stderr)
                print(res.stderr, file=sys.stderr)
                sys.exit()

        if pe:
            pardre_cmd = PARDRE.format(
                pardre_fwd_infile, pardre_rev_infile,
                dedupfile_fwd, dedupfile_rev,
                PARDRE_L, PARDRE_C,
                outprefix, outprefix,
            )
        else:
            pardre_cmd = PARDRE.format(
                pardre_fwd_infile,
                dedupfile_fwd,
                PARDRE_L, PARDRE_C,
                outprefix, outprefix,
            )
        print(f"\nRunning ParDRe command:\n{pardre_cmd}\n")
        par_res = subprocess.run(pardre_cmd, shell=True, capture_output=True)
        if par_res.returncode != 0:
            print(f"{pardre_cmd} returned non-zero exit status:", file=sys.stderr)
            print(par_res.stderr, file=sys.stderr)
            sys.exit()

        if pe:
            umi_cutadapt_cmd = f"cutadapt --quality-base={PHRED_BASE} "\
                f"-u {cut_u} -U {cut_U} --match-read-wildcards "\
                f"-o {trim_in_fwd} -p {trim_in_rev} {dedupfile_fwd} {dedupfile_rev} "\
                f"> {outprefix}_cutadapt_umi.log 2> {outprefix}_cutadapt_umi.err"
        else:
            umi_cutadapt_cmd = f"cutadapt --quality-base={PHRED_BASE} "\
                f"-u {cut_u} --match-read-wildcards "\
                f"-o {trim_in_fwd} {dedupfile_fwd} "\
                f"> {outprefix}_cutadapt_umi.log 2> {outprefix}_cutadapt_umi.err"

        print(f"\nRunning cutadapt command to remove umi:\n{umi_cutadapt_cmd}\n")
        subprocess.call(umi_cutadapt_cmd, shell=True)

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
    adapts = samp_info["adapter_seqs"]
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
            "adapseq": adapts[i],
            "phredbase": conf_dict_global["general"]["phredbase"],
            "outprefix": rep_names[i],
            "paired_reads": PE,
        }
        preprocess_file(samp_dict)
    
    # when finished with this sample's preprocessing, move back up
    #  to condition-level directory
    os.chdir(STARTDIR)

