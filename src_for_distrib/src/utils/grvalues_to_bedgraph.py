#!/usr/bin/python

# write the values from a gr file into an output bedgraph file

import ipod_utils
import sys
import argparse

BED_FMT_STRING="%s\t%i\t%i\t%s\t%f\n" # genome_name startloc endloc sitename score

parser = argparse.ArgumentParser()

parser.add_argument('gr_file', help="gr file containing the flag locations")
parser.add_argument('output_file', help="name of output .bed file")
parser.add_argument('--chrname', help="name of chromosome to use in output", default="E_coli_genome")
parser.add_argument('--sitename', help="name to use in the site name file", default="peak")
parser.add_argument('--bedspacing', help="space (in bp) between entries in the gr file", default=5, type=int)

args=parser.parse_args()


ingr=args.gr_file
outbed=args.output_file

locs, vals = ipod_utils.read_grfile(ingr)
ostr=open(outbed, 'w')

for loc, val in zip(locs, vals):
    ostr.write(BED_FMT_STRING % ( args.chrname, loc, loc+args.bedspacing, args.sitename, val) )

ostr.close()






