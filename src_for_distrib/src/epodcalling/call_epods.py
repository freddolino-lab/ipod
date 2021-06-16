#!/usr/bin/python

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

    if invert:
        mean512_file = "{}_inverted_mean512.bedgraph".format(outprefix)
        mean256_file = "{}_inverted_mean256.bedgraph".format(outprefix)
        output_epod_file = "{}_inverted_epods.bedgraph".format(outprefix)
        output_peak_file = "{}_inverted_epods.narrowpeak".format(outprefix)
        output_epod_file_strict = "{}_inverted_epods_strict.bedgraph".format(outprefix)
        output_peak_file_strict = "{}_inverted_epods_strict.narrowpeak".format(outprefix)
    else:
        mean512_file = "{}_mean512.bedgraph".format(outprefix)
        mean256_file = "{}_mean256.bedgraph".format(outprefix)
        output_epod_file = "{}_epods.bedgraph".format(outprefix)
        output_peak_file = "{}_epods.narrowpeak".format(outprefix)
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

    # calculate rolling means and write results to bedgraph files
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
        if args.invert_scores:
            scores += -1
        starts = ctg_info.fetch_array("start")
        ends = ctg_info.fetch_array("end")

        mean_256 = pu.calc_ctg_running_mean(scores, 257, res)
        mean_512 = pu.calc_ctg_running_mean(scores, 513, res)

        # add info from running mean arrays to bedgraph objects
        for i in range(len(mean_512)):
            start = starts[i]
            end = ends[i]
            score_256 = mean_256[i]
            score_512 = mean_512[i]
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
    print("Writing to {} and {}".format(mean256_file, mean512_file))
    bedgraph_256_out.write_file(mean256_file)
    bedgraph_512_out.write_file(mean512_file)
    print("------------------------------------------------\n")

    pu.identify_epods_v3_bedgraph(
        mean512_file,
        mean256_file,
        loose_len,
        epod_loose_out,
        delta = 25,
    )
    epod_loose_out.write_file(output_epod_file)

    pu.identify_epods_v3_bedgraph(
        mean512_file,
        mean256_file,
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
