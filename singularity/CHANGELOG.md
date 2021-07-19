# Changelog

All notable changes to singularity containers will be documented
in this file.

## In progress

Here are the changes to the development code not in a singularity container.

### Added

+ added control over whether we pass --no-unal to bowtie2 during alignment step
    + You now MUST have the ["alignment"]["write_unaligned_reads_to_bam"] option
     set to either true or false

## Unreleased

## 2.3.7 - 2021-07-16

### Changed

+ Removed --no-unal flag from bowtie2 call in `run_all_alignments.py`
    + This is fine as long as the user filters unaligned reads out in the bootstrap step
