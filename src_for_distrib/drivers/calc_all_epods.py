#!/usr/bin/python

# calculate all EPODs for a set of IPOD samples that I am interested in
# we read all of the instances to look at from a table of name/config file pairs

import ipod_utils
import analyze_peaks
import sys
import argparse
import os
import os.path
import toml
import subprocess
from multiprocessing import Pool
import tempfile

# set up to parse command line arguments
parser = argparse.ArgumentParser()

#parser.add_argument('--basedir',help="root directory to find all files of interest (default: current directory)", default=os.getcwd())
parser.add_argument('--outdir',help="root directory for output files (default: current directory)", default=os.getcwd())
parser.add_argument('--numproc',help="number of processors to use (default: 12)", default=12,type=int)
parser.add_argument('--scoretype',help="Type of scores to use in EPOD calling (default: _v6rz_chipsub)", default="_v6rz_chipsub")
parser.add_argument('--invert_scores',help="look for regions of very low occupancy instead of very high. Adds _LPOD to all filenames and acts on negative occupancy", action='store_true')
parser.add_argument('config_file', help='Top level configuration file for the set of samples under study')

args=parser.parse_args()

# define a helper function to actually call epods for a single condition

def do_epod_calls(gr_file_in, outprefix, genomelength):
    """
    Do all epod calling for a given gr file, writing the results to out_prefix
    """

    print "started EPOD processing for %s" % outprefix

    med512_file = "%s_median512.gr" % outprefix
    med256_file = "%s_median256.gr" % outprefix
    output_epod_file = "%s_epods_v3.gr" % outprefix
    output_peak_file = "%s_epod_locs_v3.txt" % outprefix
    output_epod_file_strict = "%s_epods_strict_v3.gr" % outprefix
    output_peak_file_strict = "%s_epod_locs_strict_v3.txt" % outprefix

    ipod_utils.do_runningavg_opt(gr_file_in, med512_file, width=513,genomelength=genomelength)
    ipod_utils.do_runningavg_opt(gr_file_in, med256_file, width=257,genomelength=genomelength)
    ipod_utils.identify_epods_v3(med512_file, med256_file, 1024, output_epod_file,delta=25)
    ipod_utils.identify_epods_v3(med512_file, med256_file, 1024, output_epod_file_strict,delta=10)

    newpeaks = analyze_peaks.read_flaggrfile(output_epod_file)
    ostr=open(output_peak_file,'w')
    for start,end in newpeaks.get_peaks():
      ostr.write("%s %s\n" % (start,end))

    ostr.close()

    newpeaks = analyze_peaks.read_flaggrfile(output_epod_file_strict)
    ostr=open(output_peak_file_strict,'w')
    for start,end in newpeaks.get_peaks():
      ostr.write("%s %s\n" % (start,end))

    ostr.close()
    print "finished EPOD processing for %s" % outprefix

## parse the config file and obtain needed information
conf_main = toml.load(args.config_file)
BASEDIR=conf_main["general"]["basedir"]
GENOMESIZE=conf_main["genome"]["genome_length"]
CONDLIST=conf_main["general"]["condition_list"]

if os.path.exists( args.outdir ):
    if not os.path.isdir(args.outdir):
        raise(NotADirectoryError)
else:
    os.mkdir( args.outdir )
    


## now actually go through the conditions of interest and run the analysis

myprocs=Pool(args.numproc)

conf_str=open(CONDLIST)
tmp_files = []

for line in conf_str:
    dirname,conffile = line.rstrip().split()
    cond_conf = toml.load(os.path.join(BASEDIR, dirname, conffile))

    if args.invert_scores:
        gr_file_in = os.path.join(BASEDIR, dirname, cond_conf["general"]["output_path"], cond_conf["general"]["out_prefix"] + args.scoretype + ".gr")
        gr_file_tmp = tempfile.NamedTemporaryFile(suffix='.gr', delete=False)
        locs,vals=ipod_utils.read_grfile(gr_file_in)
        ipod_utils.write_grfile( locs, -1*vals, gr_file_tmp.name)
        outprefix = os.path.join( args.outdir, cond_conf["general"]["out_prefix"] + args.scoretype + "_LPOD_epods")
        myprocs.apply_async(do_epod_calls, [gr_file_tmp.name, outprefix, GENOMESIZE])
        tmp_files.append(gr_file_tmp.name)

    else:
        gr_file_in = os.path.join(BASEDIR, dirname, cond_conf["general"]["output_path"], cond_conf["general"]["out_prefix"] + args.scoretype + ".gr")
        outprefix = os.path.join( args.outdir, cond_conf["general"]["out_prefix"] + args.scoretype + "_epods")
        myprocs.apply_async(do_epod_calls, [gr_file_in, outprefix, GENOMESIZE])

conf_str.close()
myprocs.close()
myprocs.join()

for f in tmp_files:
    os.remove(f)
