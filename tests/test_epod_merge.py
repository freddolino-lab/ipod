#!/usr/bin/env python3

import unittest
import os
import sys
import pathlib

this_path = pathlib.Path(__file__).parent.absolute()
bin_path = os.path.join(this_path, '../src_for_distrib/drivers')
utils_path = os.path.join(this_path, '../src_for_distrib/src/utils')
sys.path.insert(0, bin_path)
sys.path.insert(0, utils_path)

import do_peak_and_epod_calls as pepod
import anno_tools as anno

"""
Unit tests for ipod analysis pipeline
"""

class NarrowPeakDataHandlingTestCase(unittest.TestCase):
    '''A test case to test various aspects of NarrowPeakData handling'''

    def setUp(self):

        self.fnames = ['test1.narrowpeak', 'test2.narrowpeak']

        self.this_np = anno.NarrowPeakData()
        print("Reading {}".format(self.fnames[0]))
        self.this_np.parse_narrowpeak_file(self.fnames[0])

        self.that_np = anno.NarrowPeakData()
        print("Reading {}".format(self.fnames[1]))
        self.that_np.parse_narrowpeak_file(self.fnames[1])

        self.sort_target = anno.NarrowPeakData()
        print("Reading sorted.narrowpeak")
        self.sort_target.parse_narrowpeak_file('sorted.narrowpeak')

        self.np_dict,self.contigs = pepod.load_for_epod_merge(self.fnames)

    def test_np_compare(self):
        '''Tests comparison of two narrowpeakdata objects'''

        eq_np = anno.NarrowPeakData()
        print("Reading {}".format(self.fnames[0]))
        eq_np.parse_narrowpeak_file(self.fnames[0])

        self.assertEqual(self.this_np, eq_np)
        self.assertNotEqual(self.this_np, self.that_np)

    def test_sort(self):
        '''Tests sorting of NarrowPeakData object'''

        np = anno.NarrowPeakData()
        print("Reading concat.narrowpeak")
        np.parse_narrowpeak_file('concat.narrowpeak')

        self.assertNotEqual(np, self.sort_target)

        np.sort()

        self.assertEqual(np, self.sort_target)

    def test_merge_features(self):
        '''Tests merging of two NarrowPeakData objects'''

        epod1 = anno.NarrowPeakEntry()
        epod1.chrom_name = self.sort_target[0].chrom_name
        epod1.start = self.sort_target[0].start
        epod1.end = self.sort_target[1].end
        epod1.score = ((
                (self.sort_target[0].end - self.sort_target[0].start)
                + (self.sort_target[1].end - self.sort_target[1].start)
            )
            / (epod1.end - epod1.start)
            / len(self.fnames)
        )
        epod2 = self.sort_target[2].copy()
        epod2.score = 0.5
        epod3 = self.sort_target[3].copy()
        epod3.score = 0.5

        merge_target = [epod1, epod2, epod3]

        target = anno.NarrowPeakData()
        print("Reading merged.narrowpeak")
        target.parse_narrowpeak_file('merged.narrowpeak')

        for ctg_id in self.contigs:
            ctg_np = anno.NarrowPeakData()
            pepod.add_contig_records_to_np(
                ctg_np,
                self.np_dict,
                ctg_id,
            )

        self.sort_target[0].name = 0
        self.sort_target[1].name = 1
        self.sort_target[2].name = 1
        self.sort_target[3].name = 0

        self.assertEqual(ctg_np, self.sort_target)

        grp_ints,n_reps = pepod.get_grouped_intervals(ctg_np)

        grpd_ints_target = [
            [(self.sort_target[0], 0), (self.sort_target[1], 1)],
            [(self.sort_target[2], 1)],
            [(self.sort_target[3], 0)],
        ]

        self.assertEqual(grp_ints, grpd_ints_target)

        merged_epods = pepod.build_unified_epods_narrowpeak(
            grp_ints,
            n_reps,
        )

        self.assertEqual(merged_epods, merge_target)

if __name__ == '__main__':
    unittest.main()

