#!/usr/bin/python3

# calculate the kullback-leibler divergence for discretizations of the two data files given

import argparse
import numpy as np
import scipy.stats
from recombinator.iid_bootstrap import iid_bootstrap

parser = argparse.ArgumentParser()

parser.add_argument(
    "--peak_scores", help="Text file containing a single column of peak scores"
)
parser.add_argument(
    "--nonpeak_scores", help="Text file containing a single column of non-peak scores"
)
parser.add_argument(
    "--out_file", help="File to write observed entropy, upper, and lower limits"
)
parser.add_argument(
    "--lowest_bin", help="Sets the lowest score bin for categorizing scores"
)
parser.add_argument(
    "--highest_bin", help="Sets the highest score bin for categorizing scores"
)
parser.add_argument(
    "--alpha", help="Confidence limit that will be used to set the threshold for which the maximum observed KL divergence is within that threshold's confidence interval",
    default=0.05, type=float,
)
parser.add_argument(
    "--bootstrap_num", help="The number of bootstrap samples to take",
    default=1000, type=int,
)

args = parser.parse_args()

infile_1 = args.peak_scores
infile_2 = args.nonpeak_scores
out_fname = args.out_file

LO_BIN = args.lowest_bin
HI_BIN = args.highest_bin
N_B = args.bootstrap_num
ALPHA = args.alpha

data_1 = np.loadtxt(infile_1, unpack=True)
data_2 = np.loadtxt(infile_2, unpack=True)

counts_1,bins_1 = np.histogram(data_1, bins=np.arange(LO_BIN,HI_BIN))
counts_2,bins_2 = np.histogram(data_2, bins=np.arange(LO_BIN,HI_BIN))

#pylab.figure()
#
#pylab.hist( data_1, bins=bins_1, alpha=0.5, label="in peaks")
#pylab.hist( data_2, bins=bins_2, alpha=0.5, label="out of peaks" )
#pylab.legend()
#pylab.savefig( outprefix + "_hist.png")
entropy=scipy.stats.entropy( counts_1+1, counts_2+1 )
#entropy_joint=scipy.stats.entropy( counts_1+1, counts_joint+1 )


# now do bootstrapping to calculate the confidence interval

value_distr_1 = iid_bootstrap(data_1, replications=N_B, replace=True)
value_distr_2 = iid_bootstrap(data_2, replications=N_B, replace=True)

entropy_from_bootstrap=[]

for i in range(N_B):
    counts_1, bins_1 = np.histogram( value_distr_1[i,:], bins=np.arange(LO_BIN,HI_BIN))
    counts_2, bins_2 = np.histogram( value_distr_2[i,:], bins=np.arange(LO_BIN,HI_BIN))

    entropy_boot = scipy.stats.entropy( counts_1+1, counts_2+1 )
    entropy_from_bootstrap.append( entropy_boot )

# now get the confidence interval based on percentiles
cl_lo = scipy.stats.scoreatpercentile( entropy_from_bootstrap, 100*(ALPHA/2) )
cl_hi = scipy.stats.scoreatpercentile( entropy_from_bootstrap, 100*(1- (ALPHA/2)) )

