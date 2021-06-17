# Main configuration file

Configuration files are in toml format. An example main configuration file
can be found [here](example_main_config_file.conf).

Here, we describe the configurable options in the main configuration file
in the top-level directory for a given IPOD project.

We assume you are using our singularity container as described [here][sing-link]
in our discussion of the configurable options below.

## Configuration file organization

The main configuration file comprizes ten sections, each with several options
controlling its respective portion of the pipeline. They are as follow:

1. [general](#general)
    + [condition_list](#condition-list)
    + [bindir](#bindir)
    + [basedir](#basedir)
    + [rawdir](#rawdir)
    + [phredbase](#phredbase)
2. [genome](#genome)
    + [organism_name](#organism-name)
    + [organism_name_short](#organism-name-short)
    + [genome_base](#genome-base)
    + [resolution](#resolution)
3. [processing](#processing)
    + [processed_direc](#processed-direc)
    + [adapt_max_n](#adapt-max-n)
    + [trim_trailing_junk_length](#trim-trailing-junk-length)
    + [trim_sliding_window_length](#trim-sliding-window-length)
    + [trim_sliding_window_qual](#trim-sliding-window-qual)
    + [min_processed_readlen](#min-processed-readlen)
    + [r_paired_read_file_suffix](#r-paired-read-file-suffix)
    + [f_paired_read_file_suffix](#f-paired-read-file-suffix)
    + [r_unpaired_read_file_suffix](#r-unpaired-read-file-suffix)
    + [f_unpaired_read_file_suffix](#f-unpaired-read-file-suffix)
4. [alignment](#alignment)
    + [aligned_direc](#aligned-direc)
    + [min_fragment_length](#min-fragment-length)
    + [max_fragment_length](#max-fragment-length)
    + [align_threads](#align-threads)
5. [bootstrap](#bootstrap)
    + [bootstrap_direc](#bootstrap-direc)
    + [bs_suffix](#bs-suffix)
    + [output_path](#output-path)
    + [aln_retain_flags](#aln-retain-flags)
    + [aln_reject_flags](#aln-reject-flags)
    + [aln_mapq_filter](#aln-mapq-filter)
    + [bootstrap_samples](#bootstrap-samples)
    + [samtools_threads](#samtools-threads)
    + [bootstrap_threads](#bootstrap-threads)
6. [qc](#qc)
    + [qc_direc](#qc-direc)
    + [fastqc_threads](#fastqc-threads)
7. [qnorm](#qnorm)
    + [out_dset](#out-dset)
8. [quant](#quant)
    + [alpha](#alpha)
    + [debug](#debug)
    + [diagnostic_plots](#diagnostic-plots)
    + [min_percentile_chipsub_fit](#min-percentile-chipsub-fit)
    + [slope_increment_frac](#slope-increment-frac)
    + [quant_numproc](#quant-numproc)
9. [peaks](#peaks)
    + [output_path](#output-path)
    + [rz_thresholds](#rz-thresholds)
    + [log10p_thresholds](#log10p-thresholds)
    + [windowsize_bp](#windowsize-bp)
    + [nproc](#peaks-nproc)
10. [epods](#epods)
    + [output_path](#output-path)
    + [loose_epod_length](#loose-epod-length)
    + [strict_epod_length](#strict-epod-length)
    + [nproc](#epods-nproc)
11. [idr](#idr)
    + [threshold](#idr-threshold)

## General

This section of the TOML document sets options that are used for many sections
of the pipeline.

### Condition list

The `condition_list` option expects the name of the text file containing
the directory/configuration pairs for each condition. The structure of
the file is described in the main README.md file [here][cond-file-link].

### Bindir

The `bindir` option points the pipeline to the location of the python
scripts. If you're following the instructions for use of our singularity
container as written [here][sing-link], this should be set as follows:

```bash
bindir = "/src_for_distrib/src"
``` 

### Basedir

`basedir` sets the location of your top-level directory. When following the
instructions [here][sing-link], set it to:

```bash
basedir = "/ipod_data"
```

### Rawdir

If you have your fastq (or fastq.gz) files within your data directories' `raw`
directories, set this option to `"None"`.

However, often you may have too
too many sequencing samples to practially
move each fastq file to its appropriate data directory. What you can
instead do in this case is to nest a single directory containing all your
fastq (or fastq.gz) files within your *top-level directory*. We recommend simply
calling this new directory will all your data `raw`. Now, the user can
set the `rawdir` option as follows:

```bash
rawdir = "/ipod_data/raw"
```

If you have the `rawdir` option set in this way, then instead of needing
all your data files within their proper data directories, a symlink to 
`rawdir` will be created within each of your data directories, and the
pipeline will locate the appropriate raw data using that symlink.

### Phredbase

The `phredbase` option sets the quality score encoding. It should usually
be set to 33.

## Genome

Here the user can set options for where the pipeline will locate their
bowtie2 index for their reference genome and what their desired resolution
is for this experiment.

### Organism name

`organism_name` sets the organism name....

### Organism name short

`organism_name_short` sets a brief organism name.

### Genome base

`genome_base` will provide the bowtie2 index for the alignment step.
If you have used our singularity container as described [here][sing-link],
set `genome_base` as follows, substituting the common prefix of your
bowtie2 index files for \<idx_prefix\>:

```bash
genome_base = "/ref/<idx_prefix>"
```

### Resolution

Must be an integer. Adjust `resolution` to set the resolution over which
you'll be summarizing your results.

## Processing

This section of the TOML file controls arguments passed to cutadapt and
trimmomatic.

### Processed direc

Set `processed_direc` to control what directory within each data directory
will contain the trimmed reads. We usually set as follows:

```bash
processed_direc = "processed"
```

### Adapt max n

This is passed as an argument to cutadapt.
Set `adapt_max_n` to the maximum number of adapter occurrences to cut from
any single read. We usually set `adapt_max_n = 3`.

### Trim trailing junk length

Sets the minimum quality to keep a base at the trailing end of a read.
`trim_trailing_junk_length` is passed to trimmomatic within the pipeline
to provide a value for trimmomatic's `TRAILING` argument. We usually
set `trim_trailing_junk_length = 3`.

### Trim sliding window length

Sets the window size for calculating a sliding average of base quality.
The value of `trim_sliding_window_length` is passed to
trimmomatic in the SLIDINGWINDOW argument. We typically set
`trim_sliding_window_length = 4`.

### Trim sliding window qual

Sets average quality required for trimmomatic's SLIDINGWINDOW argument.
We usually use `trim_sliding_window_qual = 15`.

### Min processed readlen

Reads shorter than `min_processed_readlen` are thrown out entirely. We
usually set `min_processed_readeln = 10`.

### R paired read file suffix

The end of the filename for reverse reads that are still paired after
running trimmomatic. Something like
`r_paired_read_file_suffix = "_trim_rev_paired.fq.gz"` should work fine.

### F paired read file suffix

The end of the filename for forward reads that are still paired after
running trimmomatic. Something like
`f_paired_read_file_suffix = "_trim_fwd_paired.fq.gz"` should work fine.

### R unpaired read file suffix

The end of the filename for reverse reads that have no pair after
running trimmomatic. Something like
`r_unpaired_read_file_suffix = "_trim_rev_unpaired.fq.gz"` should work fine.

### F unpaired read file suffix

The end of the filename for forward reads that have no pair after
running trimmomatic. Something like
`f_unpaired_read_file_suffix = "_trim_fwd_unpaired.fq.gz"` should work fine.

## Alignment

Here we set the options for the alignemt stage of the pipeline.

### Aligned direc

`aligned_direc = "aligned"` will save alignments from bowtie2 to a direcotry
called "aligned" within a given data directory.

### Min fragment length

`min_fragment_length` sets the shortest possible distance between outer ends
of reads for a paired alignment. We typically use `min_fragment_length = 0`

### Max fragment length

`max_fragment_length` sets the longest possible distance between outer ends
of reads for a paired alignment. We typically use `max_fragment_length = 2000`.

### Align threads

`align_threads` sets the value of the `-p` argument to bowtie2.
The value you use will depend on the system you're using to run your analysis.

## Bootstrap

Here we set the options for bootstrapping of reads for estimating technical
noise affecting coverage within each sample. In practice, you might not
use these bootstrapped coverage estimates. However, if any sample from a
replicate is missing and you have a paired chip/ipod/input samples,
then imputation of the missing sample's coverage is performed. That imputation
is informed by the bootstrapping that was done at this step. In addition,
this step also simply counts the actual observed coverage at each position
of the genome, at the resolution you defined using the [resolution option](#resolution).

### Bootstrap direc

`bootstrap_direc` will be the directory into which your
bootstrapping results will be saved. The output will be an hdf5 file.

### BS suffix

`bs_suffix` sets the suffix of the output hdf5 file name.

### Output path

`output_path` within this section sets the output directory of the final
quantification of ipod and chip enrichment relative to input DNA.

### Aln retain flags

`aln_retain_flags` sets the -f argument to samtools view. We typically set
`aln_retain_flags = 3`, which keeps only paired alignments mapped in a
proper pair.

### Aln reject flags

`aln_reject_flags` sets the -F argument to samtools view. We usually use
`aln_reject_flags = 2308`. This will remove alignments which are unmapped,
not the primary alignment, or that are a supplementary alignment.

### Aln mapq filter

`aln_mapq_filter` sets the minimal mapping quality that should be retained
by samtools view. We usually use `aln_mapq_filter = 30`.

### Bootstrap samples

`bootstrap_samples` sets the number bootstrap replicates that will
be performed during read bootstrapping. We usually use 50 bootstrap reps.

### Samtools threads

`samtools_threads` sets the number of threads to use by samtools sort
(the @ argument to samtools sort). This value will depend on your system.

### Bootstrap threads

`bootstrap_threads` sets the number of threads to use in the actual
bootstrap sampling. The appropriate value here will depend on your system.

## QC

These options affect the behavior of fastqc.

### QC direc

`qc_direc` indicates the directory within your data directory
into which fastqc output will be saved.

### Fastqc threads

`fastqc_threads` sets the max number of threads fastqc will use.

## Qnorm

### Out dset

`out_dset` sets the name of the dataset that will be written to the hdf5 file
containing bootstrapping results. We typically just set `out_dset = "qnorm"`.
The dataset that will be written will contain quantile normalized coverage
values.

## Quant

Options in this section determine the behavior of many steps of ipod and chip
enrichment calculation and subtraction of RNAP chip contribution to IPOD
signal.

### Alpha

`alpha` sets the desired confidence interval for the results. We usually
stick with the conventional `alpha = 0.95` to get 95% confidence limits.

### Debug

This option currently does nothing and can even be omitted from your conf file.

### Diagnostic plots

`diagnostic_plots = true` will cause this step of the pipeline to save
plots that will be helpful in determining whether subtraction of RNAP's ChIP
contribution to IPOD singnals is necessary or appropriate. We recommend you
always run our pipeline with the option set to `true`. If the resulting
plots indicate that chip subtraction is not desireable, you can control
which samples get the chip contribution subtracted from them in
your [condition-level configuration files](conf-doc).

### Min percentile chipsub fit

`min_percentile_chipsub_fit` must be an integer, and is used to set the
percentile of RNAP ChIP data above which
data will be considered when inferring the contribution of
RANP ChIP to IPOD signals. We usually set it as follows:

```bash
min_percentile_chipsub_fit = 98
```

The choice of value is ultimately arbitrary here, but keep in mind that
higher values will tend to ensure that you are really looking at data
where RNAP is very abundant, and is thus likely to have a clear relationship
with IPOD signal.

### Slope increment frac

During the process of inferring the contribution of RNAP occupancy to IPOD
signals, we fit a line to infer the relationship between IPOD signal and
RNAP ChIP signal in the data displaying RNAP ChIP signal in the
`min_percentile_chipsub_fit` percentile, as described [above](#min-percentile-chipsub-fit).

*After* that initial fit, we incrementally adjust the slope of the initial fit
upward until the fraction of the data specified by `slope_increment_frac` is
below the line. As with the `min_percentile_chipsub_fit` option, the value
you choose here is ultimately arbitrary. Keep in mind that higher values
will cause this step to perform a more aggressive chip subtraction, and
lower values will subtract a lower number from the IPOD signals. The former
will bias your analysis toward false-negative protein occupancy domains,
and the latter will bias you toward false-positive protein occupancy domains
that could simply be a reflection of RNAP occupying the sites you end up
calling protein occupancy domains.

### Quant numproc

`quant_numproc` sets the number of chip subtraction processes that can be
run in parallel. The appropriate value here will depend on your system.

## Peaks

This section contains options that affect the performance and output
of the peak calling step.

### Output path

`output_path` simply sets the name of the directory into which your
peak calling results will be written. If the directory does not
exist at calltime, it will be created for you.

### rz thresholds

`rz_threshold` must be a list of robust z-score thresholds. Each threshold
will be used to call peaks above that threshold and save a corresponding
narrowpeaks output file.

### log10p thresholds

`log10p_thresholds` has the same function as rz_thresholds, but it takes
a list of log10p-values.

### Windowsize bp

`windowsize_bp` sets the rolling mean window width (in base pairs) over
which signals will be averaged during peak calling. We usually use
`windowsize_bp = 75`.

### Peaks nproc

`nproc` sets the max number of peak calling processes that can
be run in parallel.

## epods

These options affect the behavior of the EPOD calling phase of the pipeline.

### Output path

`output_path` just denotes the name of the directory into which your EPOD
calling results will be written. If the directory does not exist at
calltime, it will be created for you.

### Loose epod length

`loose_epod_length` sets the minimum length in base pairs that must be
satisfied for a location to be called a loose epod.

### Strict epod length

`strict_epod_length` sets the minimum length in base pairs that must be
satisfied for a location to be called a strict epod.

### nproc

`nproc` sets the max number of epod calling processes that can
be run in paralell.

## IDR

This section defines options that control the behavior of the IDR thresholding
for generating a set of consensus peaks or EPODS.

### IDR threshold

`threshold` sets the upper limit irreproducible discovery rate to consider
a given peak or EPOD as reproducible. We typically set this to 0.05.

[cond-file-link]: ../README.md#top-level-conditions-file
[sing-link]: ../README.md#singularity-use
[conf-doc]: condition_config.md#quant
