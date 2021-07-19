import toml
import sys

CFG_FMT = {
    'general': [
        'condition_list',
        'bindir',
        'basedir',
        'rawdir',
        'phredbase',
    ],
    'genome': [
        'organism_name',
        'organism_name_short',
        'genome_base',
        'resolution',
        'origin_location',
    ],
    'processing' : [
        'processed_direc',
        'adapt_max_n',
        'trim_trailing_junk_length',
        'trim_sliding_window_length',
        'trim_sliding_window_qual',
        'min_processed_readlen',
        'r_paired_read_file_suffix',
        'f_paired_read_file_suffix',
        'r_unpaired_read_file_suffix',
        'f_unpaired_read_file_suffix',
    ]
    'alignment' : [
        'aligned_direc',
        'min_fragment_length',
        'max_fragment_length',
        'write_unaligned_reads_to_bam',
        'align_threads',
    ],
    'bootstrap' : [
        'bootstrap_direc',
        'bs_suffix',
        'output_path',
        'aln_retain_flags',
        'aln_reject_flags',
        'aln_mapq_filter',
        'bootstrap_samples',
        'samtools_threads',
        'bootstrap_threads',
    ],
    'qc' : [
        'qc_direc',
        'fastqc_threads',
    ],
    'qnorm' : [
        'out_dset'
    ],
    'quant' : [
        'alpha',
        'diagnostic_plots',
        'min_percentile_chipsub_fit',
        'slope_increment_frac',
        'quant_numproc',
    ],
    'peaks' : [
        'output_path',
        'rz_thresholds',
        'log10p_thresholds',
        'windowsize_bp',
        'nproc',
    ],
    'epods' : [
        'output_path',
        'loose_epod_length',
        'strict_epod_length',
        'nproc',
    ],
    'idr' : [
        'threshold',
    ]
}

missing_sections = []
missing_subsections = []
cfg = toml.load(sys.argv[1])

for k,v in CFG_FMT.items():
    try:
        this_section = cfg[k]
    except KeyError:
        missing.append(k)
        continue
    for subcfg in v:
        if not subcfg in this_section:
            missing_subsections.append("[{}][{}]".format(k,subcfg))

print("Missing sections")
print(missing_sections)
print("Missing subsections")
print(missing_subsections)
