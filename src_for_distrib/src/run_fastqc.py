#!/usr/bin/python

# script to run fastqc on each input file in my sequencing manifest
# these should be inspected to make sure everything looks ok, likely before doing ANYTHING else

import subprocess
import shlex
import tempfile
import os
from multiprocessing import Pool
import argparse

# useful constants
READDIR="raw" ;# directory containing the raw reads and fastqc output

## Set up command line parsing
parser = argparse.ArgumentParser(description="Run fastqc on all files in the condition of iterest.")

parser.add_argument('--adapseq', default="AGATCGGAAGAGC")
parser.add_argument('--phredbase',default=33,type=int)
parser.add_argument('--prefix',default="")
parser.add_argument('--suffix1',default="_R1_001.fastq.gz")
parser.add_argument('--suffix2',default="_R2_001.fastq.gz")
parser.add_argument('--threads',default=12,type=int)

# parse the command line and set all values accordingly

myargs = parser.parse_args()
PHRED_BASE = myargs.phredbase
ADAP_SEQ = myargs.adapseq
COMMON_PREFIX = myargs.prefix
COMMON_SUFFIX_R1 = myargs.suffix1
COMMON_SUFFIX_R2 = myargs.suffix2
NPROC = myargs.threads


def run_fastqc_fqfiles(infile_f, infile_r):
  cmdline = "fastqc %s %s --noextract -t 12"  % (infile_f, infile_r)
  subprocess.call(cmdline,shell=True)

intab = open('seq_manifest.txt')
in_prefixes = []
for line in intab:
  in_prefix,out_prefix = line.rstrip().split()
  in_prefixes.append(in_prefix)

for inprefix in in_prefixes:
  infile_fwd = os.path.join(READDIR,COMMON_PREFIX + inprefix + COMMON_SUFFIX_R1)
  infile_rev = os.path.join(READDIR,COMMON_PREFIX + inprefix + COMMON_SUFFIX_R2)
#  try:
#    test1 = os.stat(infile_fwd)
#    test1 = os.stat(infile_rev) 
#  except:
#    print inprefix
#    print infile_fwd
#    print infile_rev
#    print 'found an error'
#  print '---'
  run_fastqc_fqfiles(infile_fwd,infile_rev)
