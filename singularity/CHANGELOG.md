# Changelog

All notable changes to singularity containers will be documented
in this file.

## Unreleased

## 2.3.7 - 2021-07-16

### Changed

+ Removed --no-unal flag from bowtie2 call in `run_all_alignments.py`
    + This is fine as long as the user filters unaligned reads out in the bootstrap step
