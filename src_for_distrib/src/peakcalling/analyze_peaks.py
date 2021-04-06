#!/usr/bin/python

# Tools for working with peak definitions

import numpy as np
from sys import argv,stdout
from bisect import bisect
import random
import ipod_utils
import tempfile
import subprocess
import copy
import argparse

class PeakList:
    """
    Class for storing information on a set of peak locations

    Really all this means is that one has a collection of ranges in the genome, possibly with associated scalar data
    This is useful for looking at genes, ipod/chip peaks, etc.

    The primary focus of this implementation is for ipod peaks

    Note that we use closed intervals on both sides and 0 indexing for the internal representation
    """

    def __init__(self):
        self.peakstarts = []
        self.peakends = []
        self.peakvals = []
        self.numpeaks = 0

    def add_peak(self, start, end, value):

        self.peakstarts.append(start)
        self.peakends.append(end)
        self.peakvals.append(value)
        self.numpeaks += 1

    def get_peaks(self,withdat=False):

        if withdat:
            return zip(self.peakstarts, self.peakends, self.peakvals)
        else:
            return zip(self.peakstarts, self.peakends)

    def diff_peaks(self,other):
        """
        Return a new PeakList object containing  peaks in self that are not in other

        Order is maintained, and sorting is not assumed
        """

        newlist = PeakList()

        for (start,end,val) in self.get_peaks(True):
            is_unique = True
            for (ostart, oend) in other.get_peaks():
                if (check_peak_overlap(start,end,ostart,oend)):
                    is_unique = False
                    break

            if (is_unique):
                newlist.add_peak(start,end,val)

        return newlist

def read_flaggrfile(infile):
    """
    Read a gr file containing a list of 0 (no peak) or 1 (peak) flags and convert to a PeakList (returned)
    """

    instr = open(infile)
    
    is_peak = False
    currstart = 0
    currend = 0

    allpeaks = PeakList()

    prev_offset = 0

    for line in instr:
        linearr = line.rstrip().split()
        offset = int( round(float(linearr[0])) + 0.2)
        flag = int( round(float(linearr[1])) + 0.2)

        if is_peak:

            if flag > 0:
                prev_offset = offset
                continue

            else:
                currend = prev_offset
                allpeaks.add_peak(currstart, prev_offset, 0)
                prev_offset = offset
                is_peak = False
                continue

        else:

            if flag > 0:
                currstart = offset
                is_peak = True
                prev_offset = offset
                continue

            else:
                prev_offset = offset
                continue

        raise Error("We should never get here")

    # Handle any peaks wrapping around the origin
    if is_peak:
        print "Warning: Truncating peak at end of file"
        allpeaks.add_peak(currstart, prev_offset, 0)

    instr.close()
    print allpeaks.numpeaks
    return allpeaks


def read_peakranges(infile):
    """
    Read a text file where each line is a pair of peak ranges

    Returns a PeakList object with all scalar values set to 0
    """

    allpeaks = PeakList()

    instr = open(infile)

    for line in instr:
        linearr = line.rstrip().split()
        start,end = [int(i) for i in linearr]
        allpeaks.add_peak(start, end, 0)

    instr.close()
    return allpeaks

def read_peaks_from_gff(infile):
    """
    Read a .gff file and return a PeakList showing the locations of features

    In the process we strip all other information
    We assume a tab-delimited gff with the start and end of each peak in the
    fourth and fifth fields

    remember that gff files are 1-indexed and have closed intervals on both ends
    """

    allpeaks = PeakList()

    instr = open(infile)

    for line in instr:
        if line[0] == '#':
                continue

        linearr = line.rstrip().split("\t")
        start = int(linearr[3]) - 1
        end = int(linearr[4]) - 1
        allpeaks.add_peak(start, end, 0)

    instr.close()
    return allpeaks

def read_peaks_from_bed(infile):
    """
    Read a .bed file and return a PeakList showing the locations of features

    In the process we strip all other information
    We assume a tab-delimited bed with the start and end of each peak in the
    second and third fields of the file

    remember that bed files are 0-indexed and have closed intervals on the left end but not the right end
    """

    allpeaks = PeakList()

    instr = open(infile)

    for line in instr:
        if line[0] == '#':
                continue

        linearr = line.split("\t")
        start = int(linearr[1])
        end = int(linearr[2]) - 1
        allpeaks.add_peak(start, end, 0)

    instr.close()
    return allpeaks

def peaks_to_array(peaks, locs):
    """
    Convert a PeakList into a set of 0/1 flags at each location in locs

    We assume that locs is sorted; peaks does not need to be
    """

    peak_arr = np.zeros_like(locs)

    for start,end in peaks.get_peaks():
        i=bisect(locs, start) - 2
        
        while locs[i] < start:

            if i >= len(locs) - 1:
                print "WARNING: Couldn't find peak %i-%i" % (start,end)
                break

            i += 1

        if i > (len(locs) - 1):
            print "WARNING: Couldn't find peak %i-%i" % (start,end)
            break

        peak_arr[i] = 1
        while (  (i < len(locs)) and (locs[i] < end) ):
            peak_arr[i] = 1
            i += 1

    return peak_arr

def find_overlapping_peaks(peaks1, peaks2):
    # Return PeakLists containing only peaks from each set that overlap the other

    locs1 = peaks1.get_peaks()
    locs2 = peaks2.get_peaks()

    overlaps1 = set()
    overlaps2 = set()

    for start1,end1 in locs1:
        for start2,end2 in locs2:

            if ( (start1 >= start2 and start1 <= end2) or (end1 >= start2 and end1 <= end2) or (start2 >= start1 and start2 <= end1) or (end2 >= start1 and end2 <= end1)):
                overlaps1.add( (start1, end1) )
                overlaps2.add( (start2, end2) )

    list1 = list(overlaps1)
    list2 = list(overlaps2)

    list1.sort(cmp = lambda a,b: cmp(a[0], b[0]))
    list2.sort(cmp = lambda a,b: cmp(a[0], b[0]))

    overlap_list1 = PeakList()
    overlap_list2 = PeakList()

    for start,end in list1:
        overlap_list1.add_peak(start, end, 0)

    for start,end in list2:
        overlap_list2.add_peak(start, end, 0)

    return(overlap_list1, overlap_list2)



def get_overlap_stats(peaks, targets, genome_length, step=5, print_output=True):
    """
    Print statistics on the overlaps of two sets of peaks

    This includes the fraction of the genome covered by each set,
        and how many of the peaks overlap
    """

    outfmt = "%30s || %20s | %20s"
    outfmt_i    = "%30s || %20i | %20i"
    outfmt_f = "%30s || %20f | %20f"

    tmpf = tempfile.NamedTemporaryFile(suffix='.txt')

    sample_peaks = peaks.get_peaks()
    target_peaks = targets.get_peaks()

    flag_targets = range(0,genome_length,step)
    sample_peak_arr = peaks_to_array(peaks, flag_targets)

    outmat = np.array( (flag_targets, sample_peak_arr) )
    np.savetxt(tmpf.name, outmat.transpose(), fmt="%i")

    target_peak_arr = peaks_to_array(targets, flag_targets)
    sample_peak_numprobes = np.sum(sample_peak_arr)
    target_peak_numprobes = np.sum(target_peak_arr)

    # Basic information
    if (print_output):
        print outfmt % ("", "Sample Peaks", "Target Peaks")
        print "-"*77
        print outfmt_i % ("# of peaks", len(sample_peaks), len(target_peaks))
        print outfmt_f % ("Fractional peak coverage", float(np.sum(sample_peak_arr)) / len(flag_targets), float(np.sum(target_peak_arr)) / len(flag_targets))

    # Overlaps

    # This is the number of *probes* (defined by the step) that overlap
    num_overlapping_probes = np.sum(np.logical_and(sample_peak_arr, target_peak_arr))
    if (print_output):
        print outfmt_f % ("Fraction overlapping probes", float(num_overlapping_probes) / sample_peak_numprobes, float(num_overlapping_probes) / target_peak_numprobes)

    # also calculate the enrichment of overlap
    if (print_output):
        print outfmt_f % ("Overlap enrichment", (float(num_overlapping_probes) / sample_peak_numprobes) / (float(np.sum(target_peak_arr)) / len(flag_targets)), (float(num_overlapping_probes) / target_peak_numprobes)/(float(np.sum(sample_peak_arr)) / len(flag_targets)) )

    # This is the number of *peaks* that at least partly overlap a peak from
    #       the other set

    overlap_sample, overlap_target = find_overlapping_peaks(peaks, targets)

    op_sample = overlap_sample.get_peaks()
    op_target = overlap_target.get_peaks()
    if (print_output):
        print outfmt % ("Overlapping peaks", "%5i (%.4f)" % (len(op_sample), float(len(op_sample)) / len(sample_peaks)), "%5i (%.4f)" % (len(op_target), float(len(op_target)) / len(target_peaks)))
        print "Jaccard distance: %f" % ( 1 - (float(num_overlapping_probes) / ( sample_peak_numprobes + target_peak_numprobes - float(num_overlapping_probes))) )

    return { "peaks_1" : len(sample_peaks), "peaks_2" : len(target_peaks), "frac_1" : float(np.sum(sample_peak_arr)) / len(flag_targets), "frac_2" : float(np.sum(target_peak_arr)) / len(flag_targets), "overlap_frac_1" : float(num_overlapping_probes) / sample_peak_numprobes, "overlap_frac_2" : float(num_overlapping_probes) / target_peak_numprobes, "overlap_count1" : len(op_sample), "overlap_count2" : len(op_target) }

def resample_tf_locs_rot(sample_peaks, target_peaks, genomelength, numsamplings=1000, exclude=1000):
        # generate random shufflings of the target peaks to test significance of overlap
        # here we do this by circular permutations of the target locations
        # exclude is a range of shuffle sizes that is not allowed

        overlap1, overlap2 = find_overlapping_peaks(sample_peaks, target_peaks)
        print "In original data, %i of %i sample peaks overlap with targets, and %i of %i target peaks overlap with sample" % (overlap1.numpeaks, sample_peaks.numpeaks, overlap2.numpeaks, target_peaks.numpeaks)

        s_overlaps = []
        t_overlaps = []

        tmp_peakfile=tempfile.NamedTemporaryFile(suffix='.gff')
        for i in range(numsamplings):

                if (i % 100 == 0):
                        print "Working on sample %i" % i

                shuff_size = random.randint(exclude,genomelength-exclude)
                target_peaks_shuffled = PeakList()

                for startloc,endloc,value in target_peaks.get_peaks(withdat=True):
                        # note that this won't work quite right for a shuffled peak that ends up straddling the origin
                        target_peaks_shuffled.add_peak( (startloc + shuff_size) % genomelength, (endloc + shuff_size) % genomelength, value)



                #print "bedtools shuffle -i %s -g %s > %s" % (target_peak_file, genomefile, tmp_peakfile.name)
                s_overlap_d, t_overlap_d = find_overlapping_peaks(sample_peaks, target_peaks_shuffled)
                s_overlaps.append(s_overlap_d.numpeaks)
                t_overlaps.append(t_overlap_d.numpeaks)
        
        print "Median values: %i sample peaks overlap with target, %i target peaks overlap with sample" % (np.median(s_overlaps), np.median(t_overlaps))
        print "Of %i samples, %i had more (or equal) sample overlap and %i had more  (or equal) target overlap than the input distribution" % (numsamplings, np.sum(np.greater_equal(s_overlaps, overlap1.numpeaks)), np.sum(np.greater_equal(t_overlaps, overlap2.numpeaks)))

def resample_tf_locs(sample_peaks, target_peak_file, genomefile, numsamplings=100):
        # generate random shufflings of the target peaks to test significance of overlap
        # the target peak file MUST be a gff

        target_peaks = read_peaks_from_gff(target_peak_file)

        overlap1, overlap2 = find_overlapping_peaks(sample_peaks, target_peaks)
        print "In original data, %i of %i sample peaks overlap with targets, and %i of %i target peaks overlap with sample" % (overlap1.numpeaks, sample_peaks.numpeaks, overlap2.numpeaks, target_peaks.numpeaks)

        s_overlaps = []
        t_overlaps = []

        tmp_peakfile=tempfile.NamedTemporaryFile(suffix='.gff')
        for i in range(numsamplings):
                if (i % 10 == 0):
                        print "Working on sample %i" % i

                #print "bedtools shuffle -i %s -g %s > %s" % (target_peak_file, genomefile, tmp_peakfile.name)
                subprocess.call("bedtools shuffle -i %s -g %s > %s" % (target_peak_file, genomefile, tmp_peakfile.name), shell=True)
                decoy_targets = read_peaks_from_gff(tmp_peakfile.name)
                s_overlap_d, t_overlap_d = find_overlapping_peaks(sample_peaks, decoy_targets)
                s_overlaps.append(s_overlap_d.numpeaks)
                t_overlaps.append(t_overlap_d.numpeaks)
        
        print "Median values: %i sample peaks overlap with target, %i target peaks overlap with sample" % (np.median(s_overlaps), np.median(t_overlaps))
        print "Of %i samples, %i had more sample overlap and %i had more target overlap than the input distribution" % (numsamplings, np.sum(np.greater(s_overlaps, overlap1.numpeaks)), np.sum(np.greater(t_overlaps, overlap2.numpeaks)))


def resample_peak_locs(sample_peaks, target_peaks, numsamplings=1000,genomelength=4641652, logfile=None, grfile=None):
    """
    Find significance of peak assignments based on resampling
    """

    random.seed()

    numpeaks_target = target_peaks.numpeaks

    if grfile is not None:
        probelocs, probevals = ipod_utils.read_grfile(grfile)
    else:
        probelocs = np.arange(genomelength)

    logstr = None
    if (logfile is not None):
        logstr = open(logfile, 'w')
        logstr.write("sample_overlaps target_overlaps sample_overlap_frac target_overlap_frac\n")


    # Lists containing lengths of overlap1 and overlap2
    all_overlap1 = []
    all_overlap2 = []

    overlap1, overlap2 = find_overlapping_peaks(sample_peaks, target_peaks)
    print "In original data, %i of %i sample peaks overlap with targets, and %i of %i target peaks overlap with sample" % (overlap1.numpeaks, sample_peaks.numpeaks, overlap2.numpeaks, target_peaks.numpeaks)

    real_overlap1 = overlap1
    real_overlap2 = overlap2

    peaklengths = []

    for start,end in sample_peaks.get_peaks():
        peaklengths.append( end-start + 1 )

    target_peak_flags = peaks_to_array(target_peaks, probelocs)

    for sample in range(numsamplings):

        if (sample % 10 == 0):
            print "Working on sample %i" % sample
        random_peaks = PeakList()
        random_flags = np.zeros(genomelength)

        for peaklen in peaklengths:
            peakstart = random.randrange(genomelength)
            while (np.any(np.greater(random_flags[peakstart:(peakstart+peaklen)],0))):
                #print "Rejecting a peak %i %i" % (peakstart, peakstart+peaklen)
                peakstart = random.randrange(genomelength)

            random_peaks.add_peak( peakstart, peakstart+peaklen-1, 0)
            random_flags[peakstart:(peakstart+peaklen)] = 1
            #print (peakstart,peakstart+peaklen-1)

        overlap1, overlap2 = find_overlapping_peaks(random_peaks, target_peaks)
        num_overlapping_probes = np.sum(np.logical_and(random_flags, target_peak_flags))

        all_overlap1.append(overlap1.numpeaks)
        all_overlap2.append(overlap2.numpeaks)

        if (logstr is not None):
            logstr.write("%i %i %f %f \n" % (overlap1.numpeaks, overlap2.numpeaks, float(num_overlapping_probes) / np.sum(random_flags), float(num_overlapping_probes) / np.sum(target_peak_flags) ))


    print "Median values: %i sample peaks overlap with target, %i target peaks overlap with sample" % (np.median(all_overlap1), np.median(all_overlap2))
    print "Of %i samples, %i had more sample overlap and %i had more target overlap than the input distribution" % (numsamplings, np.sum(np.greater(all_overlap1, real_overlap1.numpeaks)), np.sum(np.greater(all_overlap2, real_overlap2.numpeaks)))

    if (logstr is not None):
        logstr.close()
    
def check_peak_overlap(start1,end1,start2,end2):
    #Return True iff the peaks overlap
    if ( (start1 >= start2 and start1 <= end2) or (end1 >= start2 and end1 <= end2) or (start2 >= start1 and start2 <= end1) or (end2 >= start1 and end2 <= end1)):
        return True

    return False


def main(peakfile1, peakfile2, genome_length, print_output=True, do_resample=False):
    """
    Read in a pair of peak locations and print overlap statistics
    """


    if peakfile1[-3:] == ".gr":
        peaks1 = read_flaggrfile(peakfile1)
    elif peakfile1[-4:] == ".gff":
        peaks1 = read_peaks_from_gff(peakfile1)
    elif peakfile1[-4:] == ".bed":
        peaks1 = read_peaks_from_bed(peakfile1)
    else:
        peaks1 = read_peakranges(peakfile1)

    if peakfile2[-4:] == ".gff":
        print "Recognized gff"
        peaks2 = read_peaks_from_gff(peakfile2)
    elif peakfile2[-4:] == ".bed":
        print "Recognized bed"
        peaks2 = read_peaks_from_bed(peakfile2)
    elif peakfile2[-3:] == ".gr":
        peaks2 = read_flaggrfile(peakfile2)
    else:
        peaks2 = read_peakranges(peakfile2)

    get_overlap_stats(peaks1, peaks2, genome_length=genome_length, print_output = print_output)
    #resample_peak_locs(peaks1, peaks2, logfile='resample.log')
    #resample_tf_locs(peaks1, peakfile2, genomefile)
    if do_resample:
        stdout.flush()
        resample_tf_locs_rot(peaks1, peaks2, genome_length)


if (__name__ == "__main__"):

    parser = argparse.ArgumentParser()

    parser.add_argument('file1', help="first file for comparison")
    parser.add_argument('file2', help="second file for comparison")
    parser.add_argument('--do_shuffle', help="if specified, also do a resampling test", action='store_true')
    parser.add_argument('--genome_length', help="length of the chromosome to consider (default 4641652)", type=int, default=4641652)

    args=parser.parse_args()

    #NOTE: use hdf5 file to get genome/contig lengths

    #file1, file2, genomefile = argv[1:4]
    #print genomefile
    #main(file1, file2, genomefile)
    print "Comparing %s with %s" % (args.file1, args.file2)
    main(args.file1, args.file2, args.genome_length, do_resample = args.do_shuffle)
    

