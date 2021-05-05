# Introduction

This source tree contains a complete system suitable for running read processing and scoring on IPOD data sets, using version 2.0.3 of the processing pipeline (currently in revision; an outgrowth of the methods described in [this paper](https://doi.org/10.1101/2020.01.29.924811)).

This repository contains code for performing IPOD data analysis. The code is based on the methods described in [this paper](https://doi.org/10.1101/2020.01.29.924811), and has been significantly modified since prior distributions of the IPOD analysis code were prepared. 

# Installation

The analysis pipeline provided here is reliant on several excellent pieces of open source software,
and in some cases requires specific versions in order to function properly.
To simplify the process of establishing a compatible environment,
we provide a [conda](https://docs.conda.io/en/latest/) environment definition in the accompanying file `conda_environment.yml`
that will provide nearly all tools necessary for running IPOD-HR;
this can be used as a checklist (and guide to specific required versions) even if not working in a conda environment.

This pipeline requires data be configured in specific locations.
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
|           +-- control_input_rep1_R1.fq
|           +-- control_input_rep1_R2.fq
|           +-- control_input_rep2_R1.fq
|           +-- control_input_rep2_R2.fq
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
|           +-- case_input_rep1_R1.fq
|           +-- case_input_rep1_R2.fq
|           +-- case_input_rep2_R1.fq
|           +-- case_input_rep2_R2.fq
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

At the top level directory,
there must be a text file containing directory/configuration pairs.
In the tree above, this file is named `conditions.txt`, and is formatted as

```
<directory1> <directory1>.conf
<directory2> <directory2>.conf
etc...
```

In the example of the tree above, "control" would be substituted for "directory1",
and "case" would be substituted for "directory2".

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
formatting details are given in the accompanying condition_config.md file.

### Data directories

Each data directory must contain:
* A directory called `raw`.
* Within the `raw` directory, the read files referenced in the condition-level configuration file for this sample type.

# Running the analysis

Assuming that the installation prerequisites described above have been met,
and that this source code distribution is present in a directory called `{SRC_LOC}`,
the entire pipeline can be run using the python program in
`{SRC_LOC}/drivers/run_all.py`,
specifying as a single command line argument the path to the `main.conf`
configuration file. For example, assuming the current directory is
`top-level-directory` from the example tree above, running

```bash
python {SRC_LOC}/drivers/run_all.py main.conf
```

would run the entire pipeline, from processing and aligning reads through
calculating protein occupancy enrichment.

# Output files

The IPOD-HR analysis pipeline will produce several intermediate files as well as
a final set of results.
Intermediate files are typically the results of individual pipeline steps
(e.g., running `bowtie`).
The final results will be written to the directory specified in
each condition-level configuration file at `general -> out_prefix`.
Within the output directory defined in `general -> out_prefix`,
the files of typical interest are:

* `OUTPREFIX_chipsub_mean.bedgraph` --
    bedgraph file containing the estimate of the mean robust Z-scores
    after ChIP subtraction.
    This is the most commonly used output in practice.
* `OUTPREFIX_rzchipsublog10p_mean.bedgraph` --
    bedgraph file containing the log10p-scaled robust Z-scores,
    effectively yielding p-values assuming a standard normal null distribution
* `OUTPREFIX_IPOD_vs_inp_lograt_mean.bedgraph` --
    bedgraph file of the IPOD/input log ratios,
    prior to subtraction of RNA polymerase occupancy
* `OUTPREFIX_CHIP_vs_inp_lograt_mean.bedgraph` --
    bedgraph file of the RNA Polymerase ChIP/input log ratios
* If your experimental design included paired sample types
    (i.e., your chip, ipod, and input data came from the same sample of a culture),
    then for each of the above files, a file will be genereated for each replicate
    with the same information.

# Example data

To provide a test case that illustrates the full set of inputs needed to apply
the IPOD-HR analysis tools in a nontrivial case,
we distribute a bundle referred to as `IPODHR_MINIMAL_TEST.tbz` which contains
the complete directory structure and raw reads needed to process an IPOD-HR experiment;
in this case, all data are taken from WT MG1655 growing in log phase in either
rich or minimal media. The full test data set can be downloaded [here](exdata).
Users will find examples for all required configuration/input files,
and can also run the complete analysis of the test data set by entering the sample
data directory and calling.

`python {SRC_PATH}/drivers/run_all.py main.conf`

Where {SRC_PATH} indicates the location of the analysis code distributed here.

Included in our sample data distribution are gold standard files for the final
results generated by running the pipeline, obtained using our development environment.
The results obtained from a test run by calling `make diff`,
which will display the magnitude of the differences between the files generated
in your run, and those present in the gold standard files.
Both the RMSD and MAX differences **should** be zero if the software environment
has been appropriately reproduced.

Note that while the data sets used in this test case are relatively small,
they still are designed to provide a non-trivial working example,
and will likely take several hours to run on a decently powerful workstation.

## Singularity Use

Download and install singularity to your system.

To use our singularity containers, email schroedj@umich.edu for access to
our Google Shared Drive. Once you are given access to the Shared Drive
you can follow [this link][singularity-link]. The version of the singularity
container with the most up-to-date code is in the "current" folder. Older
versions of the container can be found in the "archive" folder.

Download the sif file. These files are quite large for IPOD analysis (>1.5 GB).
This file contains all the necessary components to run our IPOD pipeline.
We recommend you also download the [example data](exdata) to test the singularity
container.

To test the downloaded singularity file on our provided dataset, enter the
top-level directory. Now run the singularity
container by running the following at the command line, substituting your
container's version for \<version\>:

```bash
singularity run -B $(pwd):/data -B /run/shm ipod_<version>.sif
conda activate ipod
```

In the above lines of code `$(pwd)` represents your current working directory,
so you must enter your top-level-directory for this to run as expected.
If your host operating system is an older version of Ubuntu, for python's
multiprocessing module to work properly within the container,
`-B /run/shm` must be included. On many systems `-B /run/shm` can likely be omitted.

At this point, if your singularity container ran properly, your command prompt
should look something like `(ipod) [ipod_<version>]$`, your data will be located at
`/data`, and the source code tree is located at `/src_for_distrib`.
To run the pipeline, run:

```bash
python /src_for_distrib/drivers/run_all_driver.py /data/main.conf
```

After the run is complete,
the results should be checked using the procedure described [above](#example-data). 

# TODO: Postprocessing tools

We also include in this source code distribution the python programs needed for
key postprocessing tasks from the accompamying manuscripts,
namely those used for calling individual TF binding peaks,
calling extended protein occupancy domains (EPODs), and for plotting and consensus
clustering of TFs based on their binding profiles. Documentation for these
programs is included in the postprocessing.md file in this directory. 

[exdata]: https://drive.google.com/file/d/1xLwPvE8YA_B0rv4rKdZfpqIx9dmP-b2W/view?usp=sharing
[singularity-link]: https://drive.google.com/drive/folders/1k85Ew32F2Ek36yjEVvLss4OWv_EE7rUv?usp=sharing
