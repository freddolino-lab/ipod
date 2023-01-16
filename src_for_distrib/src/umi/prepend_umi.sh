#!/bin/bash

### for prepending the final 11 bases of reads from I1 file to R1 reads
# writes new file "umi_R1.fastq.gz"

### usage
# prepend_umi.sh <idx_file.fastq.gz> <read_file.fastq.gz>

src_dir=$3
out_dir=$4
in1=$out_dir/idx.fq
in2=$out_dir/reads.fq
zcat $1 > $in1
zcat $2 > $in2

# place final 11 bases of index onto 5-prime end of read
paste $in1 $in2 | awk -f $src_dir/prepend_umi.awk - > $out_dir/umi_read.fq
# gzip output
gzip -f $out_dir/umi_read.fq
# remove temporary files
rm $in1
rm $in2

