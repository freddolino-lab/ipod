#!/usr/bin/python

import toml
import os
import sys
import numpy as np
import scipy.signal
import argparse

import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
utils_path = os.path.join(this_path, '../utils')
sys.path.insert(0, utils_path)

import anno_tools as anno
import peak_utils as pu

def get_ctg_records(bg_dat, ctg_id):
    '''Convenience function to grab from bg_dat only records from ctg_id.
    '''

    # instantiate a BEDGraphData object to hold this contig's input data
    ctg_info = anno.BEDGraphData()

    # iterate through what's in the input file
    for record in bg_dat:
        # if this record's contig is the current contig,
        #   store it in ctg_info
        if record.filter("chrom_name", ctg_id):
            ctg_info.add_entry(record)

    return ctg_info

def insert_into_bg_data(bg_dat, ctg_dat, median_arr, ctg_id):
    '''Convenience function to insert median vals as
    records into bedgraph object.
    '''

    starts = ctg_dat.fetch_array("start")
    ends = ctg_dat.fetch_array("end")

    # add info from running median arrays to bedgraph objects
    for i in range(len(median_arr)):
        start = starts[i]
        end = ends[i]
        score = median_arr[i]
        bg_dat.addline(
            ctg_id,
            start,
            end,
            score,
        )

def roll_median(bedgraph_input, bedgraph_compare, invert, res):

    bedgraph_512_out = anno.BEDGraphData()
    bedgraph_256_out = anno.BEDGraphData()

    # calculate rolling medians and write results to bedgraph files
    ctgs = bedgraph_input.ctg_names()
    for ctg_id in ctgs:

        ctg_input = get_ctg_records(bedgraph_input, ctg_id)
        ctg_compare = get_ctg_records(bedgraph_compare, ctg_id)

        # grab array of positions and array of values from ctg_info
        scores_input = ctg_input.fetch_array("score")
        scores_compare = ctg_compare.fetch_array("score")
        if invert:
            scores_input *= -1
            scores_compare *= -1

        median_256 = pu.calc_ctg_running_median(scores_input, 257, res)
        median_512 = pu.calc_ctg_running_median(scores_compare, 513, res)

        insert_into_bg_data(
            bedgraph_256_out,
            ctg_input,
            median_256,
            ctg_id,
        )
        insert_into_bg_data(
            bedgraph_512_out,
            ctg_compare,
            median_512,
            ctg_id,
        )

    return (bedgraph_256_out, bedgraph_512_out)


def do_epod_calls(bg_infile_path, bg_comparefile_path, outprefix, res,
                  invert, loose_len, strict_len):
    '''Do all epod calling for a given gedgraph file, writing results along
    the way.

    Args:
    -----
    bg_infile_path : str
        Path to the bedgraph file containing mean scores.
    bg_comparefile_path : str
        Path to the bedgraph file containing data that will be used
        to call epods.
    outprefix : str
        Characters to prepend to output file names.
    res : int
        Resolution of original data.
    invert : bool
        If False, do not invert scores. Output will be EPODS.
        If True, invert the scores such that output rows will
        represent regions of protein occupancy depletion, i.e., inverted EPODS.
    '''

    median512_file = "{}_median512.bedgraph".format(outprefix)
    median256_file = "{}_median256.bedgraph".format(outprefix)
    output_epod_file = "{}_epods_loose.bedgraph".format(outprefix)
    output_peak_file = "{}_epods_loose.narrowpeak".format(outprefix)
    output_epod_file_strict = "{}_epods_strict.bedgraph".format(outprefix)
    output_peak_file_strict = "{}_epods_strict.narrowpeak".format(outprefix)

    # read input data into bedgraphdata object
    bedgraph_compare = anno.BEDGraphData()
    bedgraph_compare.parse_bedgraph_file(bg_comparefile_path)
    bedgraph_input = anno.BEDGraphData()
    bedgraph_input.parse_bedgraph_file(bg_infile_path)

    ctgs = bedgraph_input.ctg_names()

    # instantiate bedgraphdata and narrowpeakdata objects to collect
    #  results
    epod_loose_out = anno.BEDGraphData()
    epod_loose_np_out = anno.NarrowPeakData()
    epod_strict_out = anno.BEDGraphData()
    epod_strict_np_out = anno.NarrowPeakData()

    # calculate rolling medians and write results to bedgraph files
    bedgraph_256_out,bedgraph_512_out = roll_median(
        bedgraph_input,
        bedgraph_compare,
        invert,
        res,
    )

    print("\n================================================")
    print("Writing to {} and {}".format(median256_file, median512_file))
    bedgraph_256_out.write_file(median256_file)
    bedgraph_512_out.write_file(median512_file)
    print("------------------------------------------------\n")

    pu.identify_epods_v3_bedgraph(
        median512_file,
        median256_file,
        loose_len,
        epod_loose_out,
        delta = 25,
    )
    epod_loose_out.write_file(output_epod_file)

    pu.identify_epods_v3_bedgraph(
        median512_file,
        median256_file,
        strict_len,
        epod_strict_out,
        delta = 10,
    )
    epod_strict_out.write_file(output_epod_file_strict)
   
    for ctg_id in ctgs:
        ctg_input = anno.BEDGraphData()
        ctg_loose = anno.BEDGraphData()
        ctg_strict = anno.BEDGraphData()

        for record in epod_loose_out:
            if record.filter("chrom_name", ctg_id):
                ctg_loose.add_entry(record)
        for record in epod_strict_out:
            if record.filter("chrom_name", ctg_id):
                ctg_strict.add_entry(record)
        for record in bedgraph_512_out:
            if record.filter("chrom_name", ctg_id):
                ctg_input.add_entry(record)

        starts = ctg_loose.fetch_array("start")
        ends = ctg_loose.fetch_array("end")
        flags_loose = ctg_loose.fetch_array("score")
        flags_strict = ctg_strict.fetch_array("score")
        signal = ctg_input.fetch_array("score")

        pu.get_peaks_from_binary_array(
            ctg_id,
            starts,
            ends,
            flags_loose,
            epod_loose_np_out,
            signal, # signal now comes from the 512 median scores
        )
        pu.get_peaks_from_binary_array(
            ctg_id,
            starts,
            ends,
            flags_strict,
            epod_strict_np_out,
            signal,
        )
    epod_loose_np_out.write_file(output_peak_file)
    epod_strict_np_out.write_file(output_peak_file_strict)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--main_conf",
        help = "Name of the main configuration file for this experiment.",
        required = True,
        type = str,
    )
    parser.add_argument(
        "--in_file",
        help = "Name of bedgraph file containing IPOD enrichments or chip-subtracted IPOD enrichments.",
        required = True,
        type = str,
    )
    parser.add_argument(
        "--compare_file",
        help = "Name of bedgraph file containing IPOD enrichments of chip-subtracted IPOD enrichments. These values will be used as the baseline to compare data in --in_file against for EPOD calling.",
        required = True,
        type = str,
    )
    parser.add_argument(
        "--out_prefix",
        help = "Prefix, including path, to append final elements of output file names to.",
        required = True,
        type = str,
    )
    parser.add_argument(
        "--invert_scores",
        help="look for regions of very low occupancy instead of very high. Adds _LPOD to all filenames and acts on negative occupancy",
        action="store_true",
    )
    parser.add_argument(
        "--resolution",
        help="resolution of input data. Will also be the resolution of the outputs.",
        required=True,
        type=int,
    )

    args = parser.parse_args()

    OUTPREF = args.out_prefix
    RESOLUTION = args.resolution
    INVERT = args.invert_scores
    IN_BEDGRAPH = args.in_file
    COMPARE_BEDGRAPH = args.compare_file
    conf_dict_global = toml.load(args.main_conf)
    LOOSE_LENGTH = conf_dict_global["epods"]["loose_epod_length"]
    STRICT_LENGTH = conf_dict_global["epods"]["strict_epod_length"]

    if LOOSE_LENGTH > STRICT_LENGTH:
        sys.exit("ERROR: within the epods section of your main configuration file, loose_epod_length MUST be less than or equal to strict_epod_length. You currently have set loose_epod_length to {} and strict_epod_length to {}. Please update these options and re-run epod calling.".format(LOOSE_LENGTH, STRICT_LENGTH))


    do_epod_calls(
        IN_BEDGRAPH,
        COMPARE_BEDGRAPH,
        OUTPREF,
        RESOLUTION,
        INVERT,
        LOOSE_LENGTH,
        STRICT_LENGTH,
    )
