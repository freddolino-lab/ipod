# Changelog

All notable changes to singularity containers will be documented
in this file.

## 2.8.2 - dev

### Added

Users can now set "bootstrap\_quantity" in their main configuration
file under the "bootstrap" category. Valid values are "coverage" or
"count". If the option is omitted from the file, the default, coverage,
is used, and this is the behavior that has been used since day 1.
However, if "count" is used, a simple fragment count will be performed
instead of the fragment-length-corrected coverage the "coverage" method
performs.

Using the "count" method is useful for preparing bedgraph files for
Enicherator or IPODerator.

### Changed

In `%runscript` section of singularity def file, now unsets evnironemnt
variable `_CONDA_PYTHON_SYSCONFIGDATA_NAME`.

For users with that environment variable set, there was an import
error when importing scipy.stats.

## 2.8.1

### Fixed

Fixed bug in choosing best threshold. Default value of alpha was 0.95, but
should have been 0.05, which denotes a 95% confidence limit. The bug caused
best threshold peak calling to fail.

## 2.8.0

### Changed

NOTE: If using UMIs, you MUST use version 2.8.0. There was a
bug in the UMI handling in all prior versions which has been
corrected in v2.8.0, along with switching from ParDre to 
umi\_tools for UMI handling.

## 2.7.2 - latest stable

NOTE: If using UMIs, you MUST use version 2.8.0. There was a
bug in the UMI handling in all prior versions which has been
corrected in v2.8.0, along with switching from ParDre to 
umi\_tools for UMI handling.

### Added

You can now set a `seed` option within the `bootstrap` heading
of your main configuration file. This should ONLY be done when
testing reproducibility, and it will be useful when running
on our lab's test dataset from E. coli, by setting the seed to
42. When run on the test date with `seed = 42`, running
`make diff` will result in 0 if everything is set up correctly.

In a future version the seed option will also set the seed to the PRNG used
to impute missing replicate values, and to do bootstrapping to calculate best peak threshold, but it currently doesn't.

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

## 2.7.1

NOTE: If using UMIs, you MUST use version 2.8.0. There was a
bug in the UMI handling in all prior versions which has been
corrected in v2.8.0, along with switching from ParDre to 
umi\_tools for UMI handling.

### Added

+ Support for different R1 and R2 adapter sequences.
    + Backward-compatible with the older method that assumed only one adapter sequence, same on each read.
    + See online documnetation for details on setting different R1 and R2 adapter sequences.

### Changed

+ Multithreading option for cutadapt and trimmomatic in preprocessing step
    + place `threads` option under `preprocessing` heading
    + defaults to one thread if option not provided in main config file
+ UMI incorporation is no longer a breaking change compared to v2.6.0 of the pipeline, so in 2.7.1 if no `handle_umi` option is proveded in the `preprocessing` heading, the pipeline will default to not handling UMIs.

## 2.7.0

NOTE: If using UMIs, you MUST use version 2.8.0. There was a
bug in the UMI handling in all prior versions which has been
corrected in v2.8.0, along with switching from ParDre to 
umi\_tools for UMI handling.

### BREAKING CHANGE

+ You must now have a \[quant\]\[handle\_umi\] option in your main config file. If you do not have UMIs, set it to "false", if you do have UMI's set it to "true".
    + In the case where you have UMI's an entirely separate section of the main
        config file is used to define the read de-duplication steps. See
        https://github.com/freddolino-lab/ipod/blob/main/docs/main_config.md
        for detailed documentation.

### Added

NOTE: If using UMIs, you MUST use version 2.8.0. There was a
bug in the UMI handling in all prior versions which has been
corrected in v2.8.0, along with switching from ParDre to 
umi\_tools for UMI handling.

+ Support for read de-duplication using UMIs (see github docs for details)
    + Can be done using NEB method, where UMI is final 11 bp of I1 read
    + Can also be set to user-defined length of 5-prime end of read
+ During spike-in normalization, more precise reporting of fraction of
  coverage allocated to each contig in the reference sequence. Now reports
  for each contig, not just spike-in contig and other stuff.
+ Can now use wildcards in read file names in the condition-level configuration
  files. All matched files will be cat together, then concatenated reads
  preprocessed together.

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
