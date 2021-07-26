# Changelog

All notable changes to singularity containers will be documented
in this file.

## In progress

Here are the changes to the development code in the development
singularity container, but NOT present in the current container.

### Added

+ added functionality for spike-in normalization to the pipeline
    + this change came with many adjustments to the pipeline and
    its associated configuration files. Read the docs to find out more.
+ added control over whether we pass --no-unal to bowtie2 during alignment step
    + You now MUST have the ["alignment"]["write_unaligned_reads_to_bam"] option
     set to either true or false

## Unreleased

## 2.3.7 - 2021-07-16

### Changed

+ Removed --no-unal flag from bowtie2 call in `run_all_alignments.py`
    + This is fine as long as the user filters unaligned reads out in the bootstrap step
