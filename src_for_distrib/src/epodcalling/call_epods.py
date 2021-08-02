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
conf_dict_global = toml.load(args.main_conf)
LOOSE_LENGTH = conf_dict_global["epods"]["loose_epod_length"]
STRICT_LENGTH = conf_dict_global["epods"]["strict_epod_length"]

if LOOSE_LENGTH > STRICT_LENGTH:
    sys.exit("ERROR: within the epods section of your main configuration file, loose_epod_length MUST be less than or equal to strict_epod_length. You currently have set loose_epod_length to {} and strict_epod_length to {}. Please update these options and re-run epod calling.".format(LOOSE_LENGTH, STRICT_LENGTH))

def do_epod_calls(bg_infile_path, outprefix, res, invert, loose_len, strict_len):
    '''Do all epod calling for a given gedgraph file, writing results along
    the way.

    Args:
    -----
    bg_infile_path : str
        Path to the bedgraph file containing data that will be used
        to call epods.
    outprefix : str
        Characters to prepend to output file names.
    res : int
        Resolution of original data.
    invert : bool
        If False, do not invert scores. Output will be EPODS.
        If True, invert the scores such that output rows will
        represent regions of protein occupancy depletion.
    '''

    #if invert:
    #    median512_file = "{}_inverted_median512.bedgraph".format(outprefix)
    #    median256_file = "{}_inverted_median256.bedgraph".format(outprefix)
    #    output_epod_file = "{}_inverted_epods_loose.bedgraph".format(outprefix)
    #    output_peak_file = "{}_inverted_epods_loose.narrowpeak".format(outprefix)
    #    output_epod_file_strict = "{}_inverted_epods_strict.bedgraph".format(outprefix)
    #    output_peak_file_strict = "{}_inverted_epods_strict.narrowpeak".format(outprefix)
    #else:
    median512_file = "{}_median512.bedgraph".format(outprefix)
    median256_file = "{}_median256.bedgraph".format(outprefix)
    output_epod_file = "{}_epods_loose.bedgraph".format(outprefix)
    output_peak_file = "{}_epods_loose.narrowpeak".format(outprefix)
    output_epod_file_strict = "{}_epods_strict.bedgraph".format(outprefix)
    output_peak_file_strict = "{}_epods_strict.narrowpeak".format(outprefix)

    bedgraph_input = anno.BEDGraphData()
    bedgraph_input.parse_bedgraph_file(bg_infile_path)
    bedgraph_input.cleanup()
    ctgs = bedgraph_input.ctg_names()

    bedgraph_512_out = anno.BEDGraphData()
    bedgraph_256_out = anno.BEDGraphData()
    
    epod_loose_out = anno.BEDGraphData()
    epod_loose_np_out = anno.NarrowPeakData()
    epod_strict_out = anno.BEDGraphData()
    epod_strict_np_out = anno.NarrowPeakData()

    # calculate rolling medians and write results to bedgraph files
    for ctg_id in ctgs:
        # instantiate a BEDGraphData object to hold this contig's input data
        ctg_info = anno.BEDGraphData()
        # iterate through what's in the input file
        for record in bedgraph_input:
            # if this record's contig is the current contig,
            #   store it in ctg_info
            if record.filter("chrom_name", ctg_id):
                ctg_info.add_entry(record)

        # grab array of positions and array of values from ctg_info
        scores = ctg_info.fetch_array("score")
        if invert:
            scores *= -1

        starts = ctg_info.fetch_array("start")
        ends = ctg_info.fetch_array("end")

        median_256 = pu.calc_ctg_running_median(scores, 257, res)
        median_512 = pu.calc_ctg_running_median(scores, 513, res)

        # add info from running median arrays to bedgraph objects
        for i in range(len(median_512)):
            start = starts[i]
            end = ends[i]
            score_256 = median_256[i]
            score_512 = median_512[i]
            bedgraph_512_out.addline(
                ctg_id,
                start,
                end,
                score_512,
            )
            bedgraph_256_out.addline(
                ctg_id,
                start,
                end,
                score_256,
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
        for record in bedgraph_input:
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
            signal,
            epod_loose_np_out,
        )
        pu.get_peaks_from_binary_array(
            ctg_id,
            starts,
            ends,
            flags_strict,
            signal,
            epod_strict_np_out,
        )
    epod_loose_np_out.write_file(output_peak_file)
    epod_strict_np_out.write_file(output_peak_file_strict)

do_epod_calls(
    IN_BEDGRAPH,
    OUTPREF,
    RESOLUTION,
    INVERT,
    LOOSE_LENGTH,
    STRICT_LENGTH,
)
