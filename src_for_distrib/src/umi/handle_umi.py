#!/usr/bin/python

'''Accepts a fastq files containing I1 reads, R1 reads, and R2 reads.
Extracts the region of each index 1 read corresponding to the UMI, and prepends
it to R1 read sequence. Compares the first u+5, where u is the length of the UMI,
bases of each UMI-read chimera, throwing out duplicated reads. Writes new R1
and R2 read files, now without PCR duplicates.
'''

import argparse
import os
import sys

import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, "../utils")
sys.path.insert(0, utils_path)

if __name__ == "__main__":
    # parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--I1-file',
        help = "name of fastq (can be gzipped) file containing I1 reads."
    )
    parser.add_argument(
        '--R1-file',
        help = "name of fastq (can be gzipped) file containing R1 reads."
    )
    parser.add_argument(
        '--I1-file',
        help = "name of fastq (can be gzipped) file containing R2 reads."
    )
    args = parser.parse_args()

