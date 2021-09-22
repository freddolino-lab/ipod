# Changelog

All notable changes to singularity containers will be documented
in this file.

## 2.5.2

### Changed

+ Made merged epod output files inlcude "inverted" if we were calling epods on inverted scores.

## 2.5.1

### Changed

+ finalized method for reporting merged epod scores, pvals, and qvals.

## 2.5.0

### Changed

+ No longer using idr for epod merging, instead using our own custom merging in `src_for_distrib/src/epodcalling/merge_epods.py`

## 2.4.1

### Changed

+ Added command line option to the idr package to make it default
to not writing output when the merged peak cound is less than 20.
+ switched from `scipy.signal.argrelmax` to `scipy.signal.find_peaks` for
initial guess at locations to begin epod spreading from. The reason this
was necessary is that a prior change from calculating a rolling mean to
a rolling median caused peaks to instead become plateaus.
`scipy.signal.argrelmax` does not return an index when the max is within
a plateau. `scipy.signal.find_peaks` does return an index when there's a
tie.
+ Un-commented older code to allow us to skip a potential peak site,
since `scipy.signal.find_peaks` identifies many more potential sites
than `scipy.signal.argrelmax` did.

## 2.4.0

### Added

+ we now avoid writing nan, Inf, and -Inf values to bedgraph files.
+ adjusted the way in/out files are names in inverted epod calling.
+ added functionality for spike-in normalization to the pipeline
    + this change came with many adjustments to the pipeline and
    its associated configuration files. Read the docs to find out more.
+ added control over whether we pass --no-unal to bowtie2 during alignment step
    + You now MUST have the ["alignment"]["write_unaligned_reads_to_bam"] option
     set to either true or false

## 2.3.7 - 2021-07-16

### Changed

+ Removed --no-unal flag from bowtie2 call in `run_all_alignments.py`
    + This is fine as long as the user filters unaligned reads out in the bootstrap step
