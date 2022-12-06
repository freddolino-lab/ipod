# Main configuration file

Configuration files are in toml format. An example main configuration file
can be found [here](example_main_config_file.conf).

Here, we describe the configurable options in the main configuration file
in the top_level directory for a given IPOD project.

We assume you are using our singularity container as described [here][sing-link]
in our discussion of the configurable options below.

## Configuration file organization

The main configuration file comprizes ten sections, each with several options
controlling its respective portion of the pipeline. They are as follow:

1. [general](#general)
    + [condition_list](#condition_list)
    + [bindir](#bindir)
    + [basedir](#basedir)
    + [rawdir](#rawdir)
    + [phredbase](#phredbase)
    + [paired_reads](#paired_reads)
    + [cfu_scale_fac](#cfu_scale_fac)
2. [genome](#genome)
    + [organism_name](#organism_name)
    + [organism_name_short](#organism_name_short)
    + [genome_base](#genome_base)
    + [resolution](#resolution)
    + [spike_in_name](#spike_in_name)
3. [processing](#processing)
    + [processed_direc](#processed_direc)
    + [handle_umi](#handle_umi)
    + [threads](#threads)
    + [adapt_max_n](#adapt_max_n)
    + [trim_trailing_junk_length](#trim_trailing_junk_length)
    + [trim_sliding_window_length](#trim_sliding_window_length)
    + [trim_sliding_window_qual](#trim_sliding_window_qual)
    + [min_processed_readlen](#min_processed_readlen)
    + [r_paired_read_file_suffix](#r_paired_read_file_suffix)
    + [f_paired_read_file_suffix](#f_paired_read_file_suffix)
    + [r_unpaired_read_file_suffix](#r_unpaired_read_file_suffix)
    + [f_unpaired_read_file_suffix](#f_unpaired_read_file_suffix)
4. [umi](#umi)
    + [method](#method)
    + [length](#length)
    + [pardre_l](#pardre_l)
    + [pardre_c](#pardre_c)
    + [read](#read)
5. [alignment](#alignment)
    + [aligned_direc](#aligned_direc)
    + [min_fragment_length](#min_fragment_length)
    + [max_fragment_length](#max_fragment_length)
    + [write_unaligned_reads_to_bam](#write_unaligned_reads_to_bam)
    + [align_threads](#align_threads)
6. [bootstrap](#bootstrap)
    + [bootstrap_direc](#bootstrap_direc)
    + [bs_suffix](#bs_suffix)
    + [output_path](#output_path)
    + [aln_retain_flags](#aln_retain_flags)
    + [aln_reject_flags](#aln_reject_flags)
    + [aln_mapq_filter](#aln_mapq_filter)
    + [bootstrap_samples](#bootstrap_samples)
    + [samtools_threads](#samtools_threads)
    + [bootstrap_threads](#bootstrap_threads)
7. [qc](#qc)
    + [qc_direc](#qc_direc)
    + [fastqc_threads](#fastqc_threads)
8. [norm](#norm)
    + [qnorm_dset](#qnorm_dset)
    + [spikenorm_dset](#spikenorm_dset)
    + [clip_len_bp](#clip_len_bp)
9. [quant](#quant)
    + [alpha](#alpha)
    + [diagnostic_plots](#diagnostic_plots)
    + [min_percentile_chipsub_fit](#min_percentile_chipsub_fit)
    + [slope_increment_frac](#slope_increment_frac)
    + [quant_numproc](#quant_numproc)
10. [peaks](#peaks)
    + [output_path](#output_path)
    + [rz_thresholds](#rz_thresholds)
    + [log10p_thresholds](#log10p_thresholds)
    + [windowsize_bp](#windowsize_bp)
    + [nproc](#peaks_nproc)
11. [epods](#epods)
    + [output_path](#output_path)
    + [loose_epod_length](#loose_epod_length)
    + [strict_epod_length](#strict_epod_length)
    + [nproc](#epods_nproc)
12. [idr](#idr)
    + [threshold](#threshold)

## general

This section of the TOML document sets options that are used for many sections
of the pipeline.

### condition_list

The `condition_list` option expects the name of the text file containing
the directory/configuration pairs for each condition. The structure of
the file is described in the main README.md file [here][cond-file-link].

### bindir

The `bindir` option points the pipeline to the location of the python
scripts. If you're following the instructions for use of our singularity
container as written [here][sing-link], this should be set as follows:

```bash
bindir = "/src_for_distrib/src"
``` 

### basedir

`basedir` sets the location of your top-level directory. When following the
instructions [here][sing-link], set it to:

```bash
basedir = "/ipod_data"
```

### rawdir

If you have your fastq (or fastq.gz) files within your data directories' `raw`
directories, set this option to `"None"`.

However, often you may have
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

### phredbase

The `phredbase` option sets the quality score encoding. It should usually
be set to 33.

### paired_reads

The `paired_reads` option sets whether you have done paired-end
sequencing or not. Set to true if you have, or false if you have
single-end reads.

### cfu_scale_fac

`cfu_scale_fac` sets the scaling factor you'll be applying to the
numbers you get from spike-in normalization. It's simply a way to
avoid running into problems with loss of precision when numbers get
really tiny and just multiplies your spike-in normalized values
by this number. So if you ran experiments in which 1e8 cfu were
present in your samples, you would probably want to set
`cfu_scale_fac` to something around 1e8.

## genome

Here the user can set options for where the pipeline will locate their
bowtie2 index for their reference genome and what their desired resolution
is for this experiment.

### organism_name

`organism_name` sets the organism name....

### organism_name_short

`organism_name_short` sets a brief organism name.

### genome_base

`genome_base` will provide the bowtie2 index for the alignment step.
If you have used our singularity container as described [here][sing-link],
set `genome_base` as follows, substituting the common prefix of your
bowtie2 index files for \<idx_prefix\>:

```bash
genome_base = "/ref/<idx_prefix>"
```

### resolution

Must be an integer. Adjust `resolution` to set the resolution over which
you'll be summarizing your results.

### spike_in_name

The `spike_in_name` option is set to "None" if you have not provided
a spike-in against which to normalize coverage.

If you have provided a spike-in, you must set its "chromosome" name
in you reference fasta file here. For instance, if your spike-in DNA
was named "spike_in_vitro" you would set this option in your main
configuration file as follows:

```bash
spike_in_name = "spike_in_vitro"
```

## processing

This section of the TOML file controls arguments passed to cutadapt and
trimmomatic.

### processed_direc

Set `processed_direc` to control what directory within each data directory
will contain the trimmed reads. We usually set as follows:

```bash
processed_direc = "processed"
```

### handle_umi

Set to `true` if your reads have a unique molecular identifier (UMI).
The method by which UMIs are hanlded are set in the [UMI](#umi) section below.

### threads

Set to the number of cores to use for cutadapt and trimmomatic.
Defaults to 1 if not included in the main configuration file.

### adapt_max_n

This is passed as an argument to cutadapt.
Set `adapt_max_n` to the maximum number of adapter occurrences to cut from
any single read. We usually set `adapt_max_n = 3`.

### trim_trailing_junk_length

Sets the minimum quality to keep a base at the trailing end of a read.
`trim_trailing_junk_length` is passed to trimmomatic within the pipeline
to provide a value for trimmomatic's `TRAILING` argument. We usually
set `trim_trailing_junk_length = 3`.

### trim_sliding_window_length

Sets the window size for calculating a sliding average of base quality.
The value of `trim_sliding_window_length` is passed to
trimmomatic in the SLIDINGWINDOW argument. We typically set
`trim_sliding_window_length = 4`.

### trim_sliding_window_qual

Sets average quality required for trimmomatic's SLIDINGWINDOW argument.
We usually use `trim_sliding_window_qual = 15`.

### min_processed_readlen

Reads shorter than `min_processed_readlen` are thrown out entirely. We
usually set `min_processed_readeln = 10`.

### r_paired_read_file_suffix

The end of the filename for reverse reads that are still paired after
running trimmomatic. Something like
`r_paired_read_file_suffix = "_trim_rev_paired.fq.gz"` should work fine.

### f_paired_read_file_suffix

The end of the filename for forward reads that are still paired after
running trimmomatic. Something like
`f_paired_read_file_suffix = "_trim_fwd_paired.fq.gz"` should work fine.

### r_unpaired_read_file_suffix

The end of the filename for reverse reads that have no pair after
running trimmomatic. Something like
`r_unpaired_read_file_suffix = "_trim_rev_unpaired.fq.gz"` should work fine.

### f_unpaired_read_file_suffix

The end of the filename for forward reads that have no pair after
running trimmomatic. Something like
`f_unpaired_read_file_suffix = "_trim_fwd_unpaired.fq.gz"` should work fine.

## umi

Here we set the options to adjust behavior of UMI handling steps.

### method

The `method` option sets whether to handle UMIs as placed onto reads
when using NEB's UMI primers, or whether to handle UMIs that are
on the 5<sup>$\prime$</sup> end of either read1 or read2. In the 5<sup>$\prime$</sup> case,
the read on which the UMI is found can be set using the [read](#read)
option below.

To handle NEB UMIs, set `method = "NEB"`. To handle 5<sup>$\prime$</sup> read UMIs,
set `method = "5-prime"`.

### length

Set `length` to the length of the UMI. For NEB's UMI kit, this is 11.

### pardre_l

This option sets the value of the `-l` argument used by ParDRe. See
ParDRe documentation for details.

### pardre_c

This option sets the value of the `-c` argument used by ParDRe. See
ParDRe documentation for details.

### read

Sets which read (R1 or R2) on which a UMI can be found. For NEB's UMI
kit, the UMI is actually found in the final 11 bases of the I1, read.
So when using NEB's UMIs, we prepend those final 11 bases of each I1
read to the corresponding read identified by this option. ParDRe is
then run using the chimeric UMI/sequence reads created by this
prepending of the UMI to the sequencing read. However,
if the user has selected "5-prime" for the [method](#method)
option, ParDRe simply runs on the reads.

After ParDRe runs, the UMI is clipped from the 5<sup>$\prime$</sup> end of the
appropriate read using cutadapt's `-u` or `-U` option.

## alignment

Here we set the options for the alignemt stage of the pipeline.

### aligned_direc

`aligned_direc = "aligned"` will save alignments from bowtie2 to a direcotry
called "aligned" within a given data directory.

### min_fragment_length

`min_fragment_length` sets the shortest possible distance between outer ends
of reads for a paired alignment. We typically use `min_fragment_length = 0`

### max_fragment_length

`max_fragment_length` sets the longest possible distance between outer ends
of reads for a paired alignment. We typically use `max_fragment_length = 2000`.

### write_unaligned_reads_to_bam

`write_unaligned_reads_to_bam` accepts the values of either true or false.
It controls whether unaligned reads will be written to the bowtie2 output.
We typically set this to false, such that we aren't writing unneeded information
to our alignment files, but in cases where alignment quality is poor or in
other situations where knowing what the identity of unaligned reads is,
this flag should be set to true. NOTE: in the bootstrapping section of
the configuration file, the [aln\_reject\_flags](#aln-reject-flags)
argument can be set
such that unaligned records are filtered out of the bam file during
bootstrapping.

### align_threads

`align_threads` sets the value of the `-p` argument to bowtie2.
The value you use will depend on the system you're using to run your analysis.

## bootstrap

Here we set the options for bootstrapping of reads for estimating technical
noise affecting coverage within each sample. In practice, you might not
use these bootstrapped coverage estimates. However, if any sample from a
replicate is missing and you have a paired chip/ipod/input samples,
then imputation of the missing sample's coverage is performed. That imputation
is informed by the bootstrapping that was done at this step. In addition,
this step also simply counts the actual observed coverage at each position
of the genome, at the resolution you defined using the [resolution option](#resolution).

### bootstrap_direc

`bootstrap_direc` will be the directory into which your
bootstrapping results will be saved. The output will be an hdf5 file.

### bs_suffix

`bs_suffix` sets the suffix of the output hdf5 file name.

### output_path

`output_path` within this section sets the output directory of the final
quantification of ipod and chip enrichment relative to input DNA.

### aln_retain_flags

`aln_retain_flags` sets the -f argument to samtools view. We typically set
`aln_retain_flags = 3`, which keeps only paired alignments mapped in a
proper pair.

### aln_reject_flags

`aln_reject_flags` sets the -F argument to samtools view. We usually use
`aln_reject_flags = 2308`. This will remove alignments which are unmapped,
not the primary alignment, or that are a supplementary alignment.

### aln_mapq_filter

`aln_mapq_filter` sets the minimal mapping quality that should be retained
by samtools view. We usually use `aln_mapq_filter = 30`.

### bootstrap-samples

`bootstrap_samples` sets the number bootstrap replicates that will
be performed during read bootstrapping. We usually use 50 bootstrap reps.

### samtools_threads

`samtools_threads` sets the number of threads to use by samtools sort
(the @ argument to samtools sort). This value will depend on your system.

### bootstrap_threads

`bootstrap_threads` sets the number of threads to use in the actual
bootstrap sampling. The appropriate value here will depend on your system.

## qc

These options affect the behavior of fastqc.

### qc_direc

`qc_direc` indicates the directory within your data directory
into which fastqc output will be saved.

### fastqc_threads

`fastqc_threads` sets the max number of threads fastqc will use.

## norm

### qnorm_dset

`qnorm_dset` sets the name of the dataset containing quantile normalized
coverage information. This dataset will be written to the hdf5 file
containing bootstrapping results. We typically just set `qnorm_dset = "qnorm"`.

### spikenorm_dset

`spikenorm_dset` sets the name of the dataset containing spike-in normalized
coverage information. This dataset will be written to the hdf5 file
containing bootstrapping results. We typically just set
`spikenorm_dset = "spikenorm"`.

### clip_len_bp

`clip_len_bp` sets the number of basepairs of coverage to clip from the
ends of your spikein data prior to calculating the mean coverage for
your spike-in chromosome. We find that 100 bp is usually a reasonable length to clip
from each end of the spike-in, but we recommend that you check how alignment
to your spike-in behaves, as a different clipping length may be appropriate
for other sample types and laboratories.

## quant

Options in this section determine the behavior of many steps of ipod and chip
enrichment calculation and subtraction of RNAP chip contribution to IPOD
signal.

### alpha

`alpha` sets the desired confidence interval for the results. We usually
stick with the conventional `alpha = 0.95` to get 95% confidence limits.

### diagnostic_plots

`diagnostic_plots = true` will cause this step of the pipeline to save
plots that will be helpful in determining whether subtraction of RNAP's ChIP
contribution to IPOD singnals is necessary or appropriate. We recommend you
always run our pipeline with the option set to `true`. If the resulting
plots indicate that chip subtraction is not desireable, you can control
which samples get the chip contribution subtracted from them in
your [condition-level configuration files](conf-doc).

### min_percentile_chipsub_fit

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

### slope_increment_frac

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

### quant_numproc

`quant_numproc` sets the number of chip subtraction processes that can be
run in parallel. The appropriate value here will depend on your system.

## peaks

This section contains options that affect the performance and output
of the peak calling step.

### output_path

`output_path` simply sets the name of the directory into which your
peak calling results will be written. If the directory does not
exist at calltime, it will be created for you.

### rz_thresholds

`rz_thresholds` must be a list of robust z-score thresholds. Each threshold
will be used to call peaks above that threshold and save a corresponding
narrowpeaks output file.

### log10p_thresholds

`log10p_thresholds` has the same function as rz_thresholds, but it takes
a list of log10p-values.

### windowsize_bp

`windowsize_bp` sets the rolling mean window width (in base pairs) over
which signals will be averaged during peak calling. We usually use
`windowsize_bp = 75`.

### peaks_nproc

`nproc` sets the max number of peak calling processes that can
be run in parallel.

## epods

These options affect the behavior of the EPOD calling phase of the pipeline.

### output_path

`output_path` just denotes the name of the directory into which your EPOD
calling results will be written. If the directory does not exist at
calltime, it will be created for you.

### loose_epod_length

`loose_epod_length` sets the minimum length in base pairs that must be
satisfied for a location to be called a loose epod.

### strict_epod_length

`strict_epod_length` sets the minimum length in base pairs that must be
satisfied for a location to be called a strict epod.

### nproc

`nproc` sets the max number of epod calling processes that can
be run in paralell.

## idr

This section defines options that control the behavior of the IDR thresholding
for generating a set of consensus peaks or EPODS.

### threshold

`threshold` sets the upper limit irreproducible discovery rate to consider
a given peak or EPOD as reproducible. We typically set this to 0.05.

[cond-file-link]: ../README.md#top-level-conditions-file
[sing-link]: ../README.md#singularity-use
[conf-doc]: condition_config.md#quant
