OCC_FILE="/corexfs/czieg/210624_LrpChIP/123a_liv_log_LrpKOsub.bedgraph"
OCC_FILE_FILT=`basename $OCC_FILE .bedgraph`_filtered.bedgraph
grep U00096.3 $OCC_FILE > $OCC_FILE_FILT


for cutoff in `seq 0 0.1 8` ; do
	this_peakfile="/corexfs/czieg/210624_LrpChIP/123a_liv_log_${cutoff}_peaks.bed"
	out_prefix=`basename $this_peakfile .bed`_output
	grep U00096.3 $this_peakfile > ${out_prefix}_peaks_filtered.bed
	bedtools complement -i ${out_prefix}_peaks_filtered.bed -g ecoli.genome > ${out_prefix}_complement.bed

	bedtools map -a ${out_prefix}_peaks_filtered.bed -b $OCC_FILE_FILT -o mean -c 4 > ${out_prefix}_in_peaks_occs.bed
	bedtools map -a ${out_prefix}_complement.bed -b $OCC_FILE_FILT -o mean -c 4 > ${out_prefix}_out_of_peaks_occs.bed

	awk -F "\t" '{print $5}' ${out_prefix}_in_peaks_occs.bed > in_peak_occs.txt
	awk -F "\t" '{print $4}' ${out_prefix}_out_of_peaks_occs.bed > out_peak_occs.txt

	echo -n "for cutoff $cutoff:  "
	python calc_KL_div_opt.py in_peak_occs.txt out_peak_occs.txt output_cut${cutoff} 2> /dev/null
done

