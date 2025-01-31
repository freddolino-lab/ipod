# Condition-level config files

Each condition-level directory must contain its own configuration file. They
are in the toml format.

An example condition-level configuration file can be found
[here](example_condition_config_file.conf).

Here we describe the condition-level configuration file structure.

1. [general](#general)
    + [out_prefix](#out_prefix)
    + [sample_types](#sample_types)
2. [quant](#quant)
    + [spikein_amount](#spike_in_amount)
    + [cfu](#cfu)
    + [spikenorm_samples](#spikenorm_samples)
    + [qnorm_samples](#qnorm_samples)
    + [paired](#paired)
    + [force_onesample](#force_onesample)
    + [chipsub_numerators](#chipsub_numerators)
3. [ipod](#ipod)
    + [directory](#directory)
    + [R1_raw_files](#R1_raw_files)
    + [R2_raw_files](#R2_raw_files)
    + [adapter_seqs](#adapter_seqs)
    + [R1_adapter_seqs](#R1_adapter_seqs)
    + [R2_adapter_seqs](#R2_adapter_seqs)
    + [sample_names](#sample_names)
4. [input](#input)
    + same sub-configurations as item 3 above
5. [chip](#chip)
    + same sub-configurations as item 3 above

## general

In this section of the TOML file we set several options to control
general behavior of the pipeline.

### out_prefix

`out_prefix` sets the string of characters that will be prepended
to all output files. We ususally simply set this to the condition name.
Note that you should avoid spaces in your prefix, and instead use underscores.

### sample_types

`sample_types` must be a python list. The pipeline uses this list
to set up an iterator to iterate through each element of this list, and
then to process the files of each sample type. For example, we usually have
RNAP ChIP, IPOD, and input DNA sample. Our `sample_types` option in this case
would be:

```bash
sample_types = ["chip", "inp", "ipod"]
```

In this example, the condition-level configuration file *must* have a `chip`,
an `inp`, and an `ipod` section, each section populated with its own set
of options, described further below.

## quant

Here we set options for how quantification of ChIP and IPOD enrichment
should be handled.

### spike_in_amount

The `spikein_amount` option denotes either the mass (in ng) of spikein added
to each sample in this condition, or to the number of cfu of a biological
spike-in added to each sample in this condition.

If no spike-in was performed, set to "None".

### cfu

The `cfu` option sets the total number of colony forming units' nucleic acid
represented by each replicate for this condition. For instance, if replicate 1
contained $3.0 x 10^8 cfu/mL$ and replicate 2 contained $1.0 x 10^8 cfu/mL$
and 30 mL was sampled from each, then the total cfu for each replicate
would be $30 mL x 3.0 x 10^8 cfu/mL = 90 x 10^8 cfu$ and
$30 mL x 1.0 x 10^8 cfu/mL = 30 x 10^8 cfu$, respectively. You would then
set the `cfu` option as below:

```bash
cfu = [9000000000, 3000000000]
```

This option is only used when performing spike-in normalization, and can
be set to "None" if no spike-in was provided.

### spikenorm_samples

If spike-in normalization was performed for any sample type, include it
as a string in this list. So, if you performed a spike-in for input and
ChIP samples, you would set as follows:

```bash
spikenorm_samples = ["inp", "chip"]
```

If no samples should have spike-in normalization applied, set this option
as follows:

```bash
spikenorm_samples = []
```

Note that if you have an experiment in which some sample types should
be spike-in normalized and some should be quantile normalized,
you must include your input sample name in both the `qnorm_samples` and
the `spikenorm_samples` lists. This way enrichments can be
calculated relative to spike-in normalized inputs for the spike-in samples,
and relative to quantile-normalized inputs for the qnormed samples.

### qnorm_samples

For samples that did not have spike-in normalization, quantile normalization
must be applied. To set which sample types should be quantile normalized,
use this option as follows:

```bash
qnorm_samples = ["inp", "ipod"]
```

The above setting for the `qnorm_samples` option would perform quantile
normalization for input and ipod sample types.

If no samples should have quantile normalization applied, set this option
as follows:

```bash
qnorm_samples = []
```

Note that if you have an experiment in which some sample types should
be spike-in normalized and some should be quantile normalized,
you must include your input sample name in both the `qnorm_samples` and
the `spikenorm_samples` lists. This way enrichments can be
calculated relative to spike-in normalized inputs for the spike-in samples,
and relative to quantile-normalized inputs for the qnormed samples.

### paired

`paired` is a boolean option. If your experimental design includes sampling
from a single culture for each of your IPOD, ChIP, and input samples in a
given replicate, the set `paired = true`. If, however, your IPOD, ChIP, and
input samples are not paired, set `paired = false`. This is vitally important,
as it changes the way incorporation
of replicate-to-replicate variability is performed in the
quant, peak calling, and epod calling steps.

### force_onesample

`force_onesample` is a boolean flag. It is only used in the case where
the user has provided unpaired data and for at least one sample type, there is only
only biological replicate present. We recommend setting
`force_onesample = false`, as this is the most appropriate way to
handle acutal data. Setting to true will force the pipeline to report enrichment
scores and confidence limits, despite there being only a single biological
replicate available for at least one data type in this condition. Of course
such estimates should not be trusted, but sometimes it is useful to have them
reported anyway, in which case the user can set this option to true, bearing
in mind that extreme caution must be employed when interpreting results.

### chipsub_numerators

`chipsub_numerators` must be a python list. This option sets the samples,
from the set of [sample\_types](#sample_types)
from which the inferred contribution of RNAP occupancy to the ipod score should
be subtracted from the initial IPOD scores.

For any sample type in this list, the pipeline will infer the extent to which
RNAP occupancy (from your ChIP data) has contributed to the IPOD signal
(see [here][chipsub-main-doc] for more details). After the inference is complete
the trend in RNAP ChIP association with IPOD score will be subtracted from the
original IPOD score, resulting in a ChIP-subtracted IPOD enrichment score.

To subtract RNAP occupancy from a sample called "ipod", set this option as follows:

```bash
chipsub_numerators = ["ipod"]
```

NOTE: for all samples in the `chipsub_numerators` list, separate files containing
chip-subtracted scores *and* non-chip-subtracted scores will be written.

## ipod

Here we set options to control where the pipeline will look for "ipod" data.
We skip these sections of the documentation for input and Chip data,
because they are essentially the same options, but applied to the other types
of data.

### directory

This option identifies where the pipeline will search for this sample type's data.
For example, we typically simply use `directory = "ipod"` for the IPOD data, and
this is how we would set this option in the examply directory structure in
[README.md][main-doc].

### R1_raw_files

This option must be a list of names of files containing *forward sequencing reads*.
The proper ordering of these file names is absolutely essential in order to
properly pair sequencing data with the replicates from which they derive.
The first file name should contain forward reads
from replicate 1, the second forward reads from replicate 2, etc..
This convention must then be followed for all other sample types, i.e.,
ChIP and input, in order to ensure that replicate 1 from the IPOD data is
properly paired with replicate 1 from ChIP and input data.

### R2_raw_files

Here you must provide a list of names of files containing *reverse sequencing reads*.
See the [R1\_raw\_files](#r1_raw_files) section for requirements
on the proper ordering of replicates.

### adapter_seqs

NOTE: if you have *different* adapter sequences on your R2 3<sup>$\prime$</sup> ends
than your R1 3<sup>$\prime$</sup> ends, omit this option. Instead, use the
[R1\_adapter\_seqs](#r1_adapter_seqs) and [R2\_adapter\_seqs](#r2_adapter_seqs)
options.

If your adapter sequences are the same on the 3<sup>$\prime$</sup> ends of
both your R1 and R2 reads, set the adapter sequence here as a list.
Each element in the list must be the adapter sequence to cut from
the 3<sup>$\prime$</sup> end of the R1 and R2 reads *for each replicate*.
So if you have two replicates and each replicate has the same
adapter sequence, this option would be a list of two elements,
each element being the same adapter sequence. However, if you
used different adapter sequences for the different replicates,
you could set the sequences accordingly here.

### R1_adapter_seqs

Here you must provide a list of adapter sequences to be trimmed from
the 3<sup>$\prime$</sup> ends of your R1 reads *for each replicate*. Again, refer to
the [R1\_raw\_files](#r1_raw_files) section for requirements on proper
ordering of adapter sequences to ensure that adapters are associated
with the appropriate replicates.

### R2_adapter_seqs

Same as [R1\_adapter\_seqs](#r1_adapter_seqs), but for trimming
adapters from the 3<sup>$\prime$</sup> end of R2 reads *for each replicate*.

### sample_names

This option sets what you would like the analysis pipeline to set
each replicate's sample name to be. *NOTE:* these sample names *MUST*
inlude "rep{N}", where {N} is substituted for the replicate number
to which each sample corresponds. This is necessary, because
there are points in the pipeline where the characters "rep\d+" are
searched for in file names to associate data with the appropriate
replicates. If your files do not contain "rep{N}", they will not
be found at these steps of the analysis.

For this option in IPOD data, we simply set the list to something
like `["ipod_rep1", "ipod_rep2"]`.

## input

See options for chip ipod section above

## chip

See options for chip ipod section above

[chipsub-main-doc]: main_config.md#min-percentile-chipsub-fit
[main-doc]: ../README.md#file-preparation
