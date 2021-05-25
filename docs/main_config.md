# Main configuration file

Configuration files are in toml format. 

Here, we describe the configurable options in the main configuration file
in the top-level directory for a given IPOD project.

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
    + [windowsize_bp](#windowsize-bp)
    + [nproc](#peaks-nproc)
10. [epods](#epods)
    + [output_path](#output-path)
    + [nproc](#epods-nproc)

## General

This section of the TOML document sets options that are used for many sections
of the pipeline.

### Condition list

The `condition_list` option points the pipeline to the text file containing
the directory/configuration pairs for each condition. The structure of
the file is described [here][cond-file-link].

### Bindir

### Basedir

### Rawdir

### Phredbase

## Genome

### Organism name

### Organism name short

### Genome base

### Resolution

## Processing

### Processed direc

### Adapt max n

### Trim trailing junk length

### Trim sliding window length

### Trim sliding window qual

### Min processed readlen

### R paired read file suffix

### F paired read file suffix

### R unpaired read file suffix

### F unpaired read file suffix

## Alignment

### Aligned direc

### Min fragment length

### Max fragment length

### Align threads

## Bootstrap

### Bootstrap direc

### BS suffix

### Output path

### Aln retain flags

### Aln reject flags

### Aln mapq filter

### Bootstrap samples

### Samtools threads

### Bootstrap threads

## QC

### QC direc

### Fastqc threads

## Qnorm

### Out dset

## Quant

### Alpha

### Debug

### Diagnostic plots

### Min percentile chipsub fit

### Slope increment frac

### Quant numproc

## Peaks

### Output path

### Windowsize bp

### Peaks nproc

## epods

### Output path

### epods nproc

[cond-file-link]: ../README.md#top-level-conditions-file
