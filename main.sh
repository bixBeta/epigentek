trim(){

	mkdir TrimQC_stats fastQC trimmed_fastqs
	for i in fastqs/*.gz
	do
		trim_galore --nextseq 20 --gzip --length 50  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
	done
	mv *_trimming_report.txt TrimQC_stats
	mv *trimmed.fq.gz trimmed_fastqs

}





trim
