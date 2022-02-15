# Introduction

This source tree contains code suitable for running read processing and scoring of protein occupancy in IPOD data sets using version 2.4.0 of the processing pipeline (version 1 is currently in revision; an outgrowth of the methods described in [this paper](https://doi.org/10.1101/2020.01.29.924811)).

This repository contains code for performing IPOD data analysis. The code is based on the methods described in [this paper](https://doi.org/10.1101/2020.01.29.924811), but has been significantly modified since prior distributions of the IPOD analysis code were prepared. 

# Installation

The analysis pipeline provided here is reliant on several excellent pieces of open source software,
and in some cases requires specific versions in order to function properly.
To simplify the process of establishing a compatible environment,
we provide a [conda](https://docs.conda.io/en/latest/) environment definition in the accompanying file `conda_environment.yml`
that will provide nearly all tools necessary for running IPOD-HR;
this can be used as a checklist (and guide to specific required versions) even if not working in a conda environment.

This pipeline requires that your data be configured in specific locations.
Additionally, certain file naming conventions must be followed,
and a configuration file must be present in order to run the pipeline.
For instructions on setting up the data locations and file names,
see the [File preparation](#file-preparation) section below.

To aid in the distribution and use of our pipeline, we provide a singularity container.
For instructions on use of our singularity container, see the [Singularity Use](#singularity-use) section below.

# File preparation

The IPOD processing pipeline requires a strict file hierarchy in order to properly find and process the data files.

The hierarchy should match the following scheme for an expreiment with two replicates
for each condition (rep1 and rep2):

```
top-level-directory
+-- main.conf
+-- conditions.txt
+-- control
|   +-- control.conf
|   +-- chip
|       +-- raw
|           +-- control_chip_rep1_R1.fq
|           +-- control_chip_rep1_R2.fq
|           +-- control_chip_rep2_R1.fq
|           +-- control_chip_rep2_R2.fq
|   +-- ipod
|       +-- raw
|           +-- control_ipod_rep1_R1.fq
|           +-- control_ipod_rep1_R2.fq
|           +-- control_ipod_rep2_R1.fq
|           +-- control_ipod_rep2_R2.fq
|   +-- inp
|       +-- raw
|           +-- control_input_rep1_R1.fq
|           +-- control_input_rep2_R2.fq
|           +-- control_input_rep2_R1.fq
|           +-- control_input_rep2_R2.fq
+-- case
|   +-- case.conf
|   +-- chip
|       +-- raw
|           +-- case_chip_rep1_R1.fq
|           +-- case_chip_rep1_R2.fq
|           +-- case_chip_rep2_R1.fq
|           +-- case_chip_rep2_R2.fq
|   +-- ipod
|       +-- raw
|           +-- case_ipod_rep1_R1.fq
|           +-- case_ipod_rep1_R2.fq
|           +-- case_ipod_rep2_R1.fq
|           +-- case_ipod_rep2_R2.fq
|   +-- inp
|       +-- raw
|           +-- case_input_rep1_R1.fq
|           +-- case_input_rep1_R2.fq
|           +-- case_input_rep2_R1.fq
|           +-- case_input_rep2_R2.fq
```

The hierarchy above and the listed files are further explained in the sub-sections below.

## Top level directory

### Top level configuration file

At the top level directory,
there must be a configuration file to set options for each step in the pipeline.
In the example tree above, the file is calles `main.conf`. It is a
[TOML](https://toml.io/en/) document. The user should familiarize themselves with
each of the configurable options in the top-level configuration file.
These options are documented in more detail in the
[main_config.md file][main-cfg-doc].

### Top level conditions file

Also in the top level directory, 
there must be a text file containing directory/configuration pairs.
In the tree above, this file is named `conditions.txt`, and is formatted as

```
<directory1> <directory1>.conf
<directory2> <directory2>.conf
etc...
```

In the example of the tree above, "control" would be substituted for "\<directory1\>",
and "case" would be substituted for "\<directory2\>".

All of the specified directories must exist in the top-level directory,
and must contain a configuration file corresponding to the config name given.
We refer to each of the specified subdirectories as *condition directories*
in the documentation below; typically each one will represent all of the
data present on a single biological condition.

## Condition directories

Each condition directory must contain a condition-level configuration file
(as noted above as conditions.txt), as well as three subdirectories named
`ipod`, `inp`, and `chip`.
We collectively refer to those subdirectories as the *data directories*;
they correspond to the IPOD-HR, input, and RNA polymerase ChIP-seq samples
(ipod, inp, and chip directories, respectively, in the tree above),
respectively, under the biological condition being considered.

### Condition-level configuration file

The condition-level configuration file is a [TOML](https://toml.io/en/) document;
formatting details are given in the
accompanying [condition_config.md file][cond-cfg-doc].

### Data directories

Each data directory must contain:
* A directory called `raw`.
* Within the `raw` directory, there must exist either the read files referenced in the condition-level configuration file for this sample type, or symlinks referenced in the condition-level configuration file which point to the actual read files which may reside elsewhere.

# Running the analysis

Assuming that the installattion prerequisites described above have been met,
that your current working directory is the top level directory of your project
(i.e., the directory containing the main configuration file),
and that this source code distribution is present in a directory called `{SRC_LOC}`,
the entire pipeline can be run using the python program in
`{SRC_LOC}/drivers/run_all_driver.py`,
specifying as a single command line argument the path to the `main.conf`
configuration file. For example, assuming the current directory is
`top-level-directory` from the example tree above, running

```bash
python {SRC_LOC}/drivers/run_all_driver.py main.conf
```

would run the steps of the pipeline from processing and aligning reads through
calculating protein occupancy enrichment.

Likewise, running

```bash
python {SRC_LOC}/drivers/do_peak_and_epod_calls.py main.conf
```

will run peak and EPOD calling.

Note that for each of `run_all_driver.py` and `do_peak_and_epod_calls.py`, discrete
steps of the pipeline in each script may be skipped. The help documentation
for each script contains a brief description of which steps may be skipped,
and how to skip them. The help documentation can be accessed by invoking either
script followed by `-h`. For instance:

```bash
python {SRC_LOC}/drivers/run_all_driver.py -h
```

Will display the command line options that can be passed to `run_all_driver.py` and
a brief description of each.

# Output files

The IPOD-HR analysis pipeline will produce several intermediate files as well as
a final set of results.
Intermediate files are typically the results of individual pipeline steps
(e.g., running `bowtie2`).

## Quant step output

The "quantification" step writes numerous bedgraph files containing each
containing some quantity at each position in the reference genome for
the chosen resolution for your analysis. The final results will be written
to the directory specified in each condition-level configuration file
at `general -> out_prefix`. Typical files of interest include:

* `OUTPREFIX_IPOD_vs_inp_lograt_mean.bedgraph` --
    bedgraph file containing the estimate of the mean log2-transformed
    IPOD signal enrichment over input DNA.
* `OUTPREFIX_IPOD_vs_inp_rzlograt_mean.bedgraph` --
    bedgraph file containing the estimate of the mean robust z-score for
    IPOD enrichment over input DNA.
* `OUTPREFIX_IPOD_vs_inp_rzlogratlog10p_mean.bedgraph` --
    same information as OUTPREFIX_IPOD_vs_inp_rzlograt_mean.bedgraph, but
    the values have been transformed to represent -log10(p-values).
* Confidence limits for the above estimates, calculated by jacknife sampling,
    at the confidence level specified in `quant -> alpha` in your main
    configuration file.
* If your experimental design included paired sample types,
    i.e., your chip, ipod, and input data came from the same sample of a culture,
    then for each of the above files, a file will be generated for each replicate
    with the same information.
* All the same files as above, but for your RNAP ChIP-seq data. "IPOD" in the
    file names will be substituted with CHIP for the RNAP ChIP-seq files.
* If you are subtracting the contribution of RNAP from your IPOD signal,
    similar files as above, but containing the ChIP-subtracted IPOD signal
    will be generated. These files follow the same naming convention used for
    the above files, but will contain "chipsub" in their names.

## Peak calling output

We typically call peaks at several thresholds. In the case that you have 
paired data, peaks are called in each replicate and peaks from individual
replicates are pooled
using the irreproducible discovery rate (IDR). Files are written to the
location specified by `peaks -> output_path` in the main configuration file.
In the file names below `<val>` is substituted with the threshold above
which peaks were called for that file. Note that if you did not perform
RNAP signal subtraction, "rzchipsub" will be substituted in the file names
below with "chipsub". Additionally, if you called peaks
in your RNAP ChIP-seq data, "IPOD" will be substututed with "CHIP".
Useful files output by the peak caller are:

* `OUTPREFIX_IPOD_rzchipsub_rep1_cutoff_<val>_peaks.narrowpeak` --
    The peaks identified in replicate 1. A file like this will be
    generated for each replicate.
* `OUTPREFIX_IPOD_rzchipsub_cutoff_<val>_idr_passed.narrowpeak` --
    This file contains only the peaks which had an IDR below
    the threshold denoted by `idr -> threshold` in the main configuration
    file.
* `OUTPREFIX_IPOD_rzchipsub_rep2_cutoff_<val>_peaks_vs_wt_m9_min_glu_merged_IPOD_rzchipsub_rep1_cutoff_<val>_peaks_idr.narrowpeak` --
    This file is generated by the IDR pipeline, and it contains the merged
    peaks for the indicated replicates. A file like this will be generated
    for each pair-wise comparison of replicates. Information in the file
    includes the IDR calculated for each peak. The information in this file
    is used to determine which peaks pass the IDR threshold in each pair-wise
    comparison of replicates, and that information is subsequently used
    to inform which peaks to write to the `*_idr_passed.narrowpeak` file
    above.
* `OUTPREFIX_IPOD_rzchipsub_rep2_cutoff_<val>_peaks_vs_wt_m9_min_glu_merged_IPOD_rzchipsub_rep1_cutoff_<val>_peaks_idr.narrowpeak.png` --
    This image contains four plots that can be useful in diagnosing potential
    issues that arise in IDR calculation, and are helpful in setting an
    appropriate IDR threshold.

## EPOD calling output

Files saved during EPOD calling are written to the location specified by
`epods -> output_path` in the main configuration file. Note that if you
did not subtract the RNAP contribution to your IPOD signal, "rzchipsub"
will be substituted with some other value. Also, for each file with
"rep1" in its name you should expect to see one additional file for each
replicate's EPODS. The files include:

* Paired data
    * `OUTPREFIX_IPOD_rzchipsub_rep1_epods_loose.narrowpeak` --
        Each line represents an EPOD called under the "loose" conditions
        specified in your main configuration file.
    * `OUTPREFIX_IPOD_rzchipsub_rep1_epods_strict.narrowpeak` --
        Each line represents an EPOD called under the "strict" conditions
        specified in your main configuration file.
    * `OUTPREFIX_IPOD_rzchipsub_rep1_median256.bedgraph` --
        Contains values from convolving a 256-bp-wide rolling median over
        the signal of interest. These values are used to determine an
        appropriate threshold above which to call EPODs in the data in
        the next file below.
    * `OUTPREFIX_IPOD_rzchipsub_rep1_median512.bedgraph` --
        Values arrived at by convolving a 512-bp-wide rolling median over
        the signal of interest.
    * `OUTPREFIX_IPOD_rzchipsub_loose_merged_epods.narrowpeak` --
        **This will often be a major file of interest in EPOD calling,
        as its contents will be used to filter EPODs based on how well-reproduced
        they are across replicates.**
        This file contains the results of merging EPODS from all replicates.

        * **Column 7** is the weighted mean signal for each merged EPOD. Here,
            the weighting is proportional to each replicate's EPODs' contributions
            to the length of the merged EPOD.
            Note that in a genome browser like IGB, column 7 will be called the
            "signal" for each EPOD.
        * **Column 8** is a simple fraction of replicates in which *at least some*
            portion of each merged EPOD was represented.
            Note that in most genome browsers, certainly in IGB, column 8 is
            called the "p-value". However, the numbers in column 8 are NOT p-values.
        * **Column 9** -- **NOTE**: This is usually the column of interest
            for filtering EPODs based on their reproducibility across replicates.
            Column 9 is a measure of how well-represented this EPOD was
            across all replicates. It is essentially the EPOD's mean representation.
            Consider a case of having two replicates: if an EPOD was called in only
            a single replicate, column 9 would be 0.5. In other words, for all the
            genome positions covered by this EPOD, half of the replicates supported
            the EPODs existence.
            Now consider a genomic locus in which an EPOD was called in replicate
            1 AND 2. However, for replicate 2, although the EPOD was contained
            enirely within replicate 1's EPOD call, replicate 2's EPOD call
            was exactly half the
            the length of the EPOD from replicate 1. In this case, column 9
            would be 0.75, because on average across replicates, 75% of the
            positions in this merged EPOD were called an epod in all replicates.
            Note that in most genome browsers, certainly in IGB, column 9 is
            called the "q-value". However, column 9 in this file in NOT a
            q-value.
    * `OUTPREFIX_IPOD_rzchipsub_strict_merged_epods.narrowpeak` --
        **This will often be a major file of interest in EPOD calling,
        as its contents will be used to filter EPODs based on how well-reproduced
        they are across replicates.**
        See the description for the file above.
* Both paired and unpaired data
    * `OUTPREFIX_IPOD_rzchipsub_mean_epods_loose.narrowpeak` --
        The result of EPOD calling on the estimate of the mean signal
        at the "loose" threshold.
    * `OUTPREFIX_IPOD_rzchipsub_mean_epods_strict.narrowpeak` --
        The result of EPOD calling on the estimate of the mean signal
        at the "strict" threshold.

# Example data

To provide a test case that illustrates the full set of inputs needed to apply
the IPOD-HR analysis tools in a nontrivial case,
we distribute a bundle referred to as `IPODHR_MINIMAL_TEST.tbz` which contains
the complete directory structure and raw reads needed to process an IPOD-HR experiment;
in this case, all data are taken from WT MG1655 growing in log phase in either
rich or minimal media. The full test data set can be downloaded [here][exdata].
Users will find examples of all required configuration/input files,
and can also run the complete analysis of the test data set by entering the example
data directory and calling.

`python {SRC_PATH}/drivers/run_all_driver.py main.conf`

Where {SRC_PATH} indicates the location of the analysis code distributed here.

## Testing reproducibility

NOTE: This section requires updates to reflect incorporation of peak and epod
calling, and some major updates regarding bootstrapping and spike-in normalization

Included in our sample data distribution are gold standard files for the final
results generated by running the pipeline, obtained using our development environment.
The results obtained from a test run by running `make diff`,
which will display the magnitude of the differences between the files generated
in your run, and those present in the gold standard files.
Both the RMSD and MAX differences **should** be zero if the software environment
has been appropriately reproduced.

Note that while the data sets used in this test case are relatively small,
they still are designed to provide a non-trivial working example,
and will likely take several hours to run on a decently powerful workstation.

# Postprocessing tools

We also include in this source code distribution the python programs needed for
key postprocessing tasks from the accompamying manuscripts,
namely those used for calling IPOD peaks and
calling extended protein occupancy domains (EPODs). Documentation for these
tools is included in the [postprocessing.md file][postproc-doc]. 

# Singularity Use

Download and install singularity to your system.

To use our singularity container, email schroedj@umich.edu for access to
our Google Shared Drive. Once you are given access to the Shared Drive
you can follow [this link][singularity-link]. The version of the singularity
container with the most up-to-date code is in the "current" folder. Older
versions of the container can be found in the "archive" folder.

Download the current sif file. It will be quite large for IPOD analysis (>2 GB).
This file contains all the necessary components to run our IPOD pipeline.
We recommend you also download the [example data][exdata] to test the singularity
container.

To test the downloaded singularity file on our provided dataset, enter the
top-level directory for the example data. Now run the singularity
container by running the following at the command line, substituting your
container's version for \<version\> and substituting the **absolute path**
to the directory containing your bowtie2 index for \<path/to/ref/direc\>:

```bash
cd <top-level-directory>
singularity run \
    -B $(pwd):/ipod_data \
    -B <path/to/ref/direc>:/ref:ro \
    -B /run/shm \
    -B <path/to/raw/data/direc>:<path/to/raw/data/direc>:ro \ 
    ipod_<version>.sif
```

In the above lines of code, substitute your top-level directory for
\<top-level-directory\> to move into the top-level directory for your
experiment. `$(pwd)` represents your current working directory,
so you must first enter your top-level-directory for this to run as expected.
If your host operating system is an older version of Ubuntu, for python's
multiprocessing module to work properly within the container,
`-B /run/shm` must also be included.
On many systems `-B /run/shm` can likely be omitted.
`-B <path/to/raw/data/direc>:<path/to/raw/data/direc>:ro`
is only necessary if your `raw` directories
within each sample directory are going to be symbolic links to the actual
directory containing your raw data. If you have your actual raw data
in the raw directories, as drawn in the directory tree above, this line
should not be included. Note that when you are symlinking to
\<path/to/raw/data/direc\>, the pipeline will set up the symlinks for you,
provided you have appropriately set the location to which the symlinks
point in the main configuration file's `rawdir` option (see [here][raw-doc]
for the documentation for this option).

At this point, if your singularity container ran properly, your command prompt
should look something like `(ipod) [ipod_<version>]$`, your data will be located at
`/data`, and the source code tree is located at `/src_for_distrib`.
To run the pipeline, run:

```bash
python /src_for_distrib/drivers/run_all_driver.py /ipod_data/main.conf
```

After this first major step of the pipeline has completed running,
you may want to call peaks and epods ([see here for more information regarding peak and epod calling](#peak-epod-calls)), which can be done within the
singularity conatiner if the following manner:

```bash
python /src_for_distrib/drivers/do_peak_and_epod_calls.py /ipod_data/main.conf
```

Once this step is complete,
the results should be checked using the procedure
described [above](#testing-reproducibility). 

[exdata]: https://drive.google.com/drive/folders/1wM0EL99ypczDJJn9n-Hpz1zQF_8EmIDV?usp=sharing
[singularity-link]: https://drive.google.com/drive/folders/1ZxtYSBBaKPQAxMzOF9hmf2ec7epTHaSf?usp=sharing
[main-cfg-doc]: docs/main_config.md
[cond-cfg-doc]: docs/condition_config.md
[postproc-doc]: docs/postprocessing.md
[raw-doc]: docs/main_config.md#rawdir
[peak-epod-calls]: docs/postprocessing.md
