# Changelog

All notable changes to singularity containers will be documented
in this file.

## 2.7.2

### Changed

Fixed a bug in the quant step. The bug affected the way output file
names were generated. The bug would stop at the
first instance of rep<N> when looking for replicate ids, and so
an instance of rep<N> anywhere in the path other than the basename
would end up being used as the replicate id, rather than the desired
behavior of the instance of rep<N> in the basename.

Now the regular expression search to name
the replicate ids in output file names happens at the replicate
file basenames, not over the entire path. 

## 2.7.0

### Added

+ Support for read de-duplication using UMIs (see github docs for details)
    + Can be done using NEB method, where UMI is final 11 bp of I1 read
    + Can also be set to user-defined length of 5-prime end of read

## 2.6.0

### BREAKING CHANGE

+ Switched condition-level configuration option \[quant\]\[force\_onesample\_unpaired\] to
\[quant\]\[force\_onesample\].

### Changed

+ For multi-contig references, now doing median normalization separately for each contig, NOT globally.

## 2.5.7

### Changed

+ Fixed a bug introduced in v2.5.5 by they way I was handling thresholds at which no peaks were identified. The bug caused later EPOD calling to fail due to having re-used a variable name by accident.

## 2.5.6

### Changed

+ Spike-in normalization is now automatically skipped if the "spikin_name"
option is "None". No need to supply `--skipsteps spikenorm` at the command
line.

## 2.5.5

### Added

+ Pseudocounts are added to both spikein-normalized, and quantile-normalized, zero counts. Prior to this version, they were only usef in spike-in normalization.
+ Best peak calling threshold is determined using KL-divergence between signals in peak vs. non-peak regions.

## 2.5.4

### Added

+ Peak calling now calls peaks in RNAP-chip data

### Changed

+ adjusted the way a range of the genome is searched for epods so that we no longer wrap around the ends of contigs. Although for circular chromosomes wrapping is appropriate, we were finding we had issues when wrapping around plasmids. The issue could cause us to search in regions of the plasmid that were disjoint from the region of the plasmid actually containing our epod.
+ adjusted peak calling so that when no peaks are calling in a given replicate at a given threshold, an empty file is generated to keep later IDR calculation step from producing an error that prevents later stages of the pipeline from running.

## 2.5.3

### Changed

+ Switched to bowtie2.4.4. 
+ Added manual creation of tmpdir so bowtie2 will not run into issues with permissions when it sets up FIFOs.

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
