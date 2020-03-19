trim(){

	mkdir TrimQC_stats fastQC trimmed_fastqs
	for i in fastqs/*.gz
	do
		trim_galore --quality 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
	done
	mv *_trimming_report.txt TrimQC_stats
	mv *trimmed.fq.gz trimmed_fastqs

}


DIR='/Volumes/TOSHIBA/SQ190110D–R6_Colombia/Analysis/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa'

align(){

	cd trimmed_fastqs

	for i in *.gz 
	do
		iSUB=`echo $i | cut -d "_" -f1`
		bwa mem -t 4 -M -R "@RG\tID:${iSUB}\tSM:${iSUB}\tPL:ILLUMINA\tLB:${iSUB}\tPU:1" ${DIR} $i \
                | samtools view -@ 4 -b -h -F 0x0100 -O BAM -o ${iSUB}.bam
    done

    mkdir primary-BAMS
	mv *.bam primary-BAMS
	mv primary-BAMS ..

	cd ..


}

sort(){
		cd primary-BAMS

		    for i in *.bam
		    do
		    samtools sort $i > `echo  $i | cut -d "." -f1`.sorted.bam
		    done

		    for i in *.sorted.bam
		    do
		    	samtools index $i
		    done

			# alignment stats etc. on raw bams
			for i in *.sorted.bam
			do
				iSUB=`echo $i | cut -d "." -f1`
				samtools flagstat $i > ${iSUB}.primary.flagstat
				samtools idxstats $i > ${iSUB}.primary.idxstats
			done

		cd ..

		
}



markDups(){
		

		cd primary-BAMS
        
        for i in *.sorted.bam
        do
            iSUB=`echo $i | cut -d "." -f1`
            
            picard \
            MarkDuplicates \
            INPUT=$i \
            OUTPUT=${iSUB}.dupMarked.sorted.bam \
            ASSUME_SORTED=true \
            REMOVE_DUPLICATES=false \
            METRICS_FILE=${iSUB}.MarkDuplicates.metrics.txt \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR=tmp

        done
        mkdir dupMarked-BAMS
        mv *.dupMarked.sorted.bam dupMarked-BAMS

		cd ..

}




# alignBowtie2(){

# 	cd trimmed_fastqs

# 	for i in *.gz 
# 	do
# 		iSUB=`echo $i | cut -d "_" -f1`
		
# 		(bowtie2 \
# 		--local \
#         -x /Volumes/TOSHIBA/SQ190110D–R6_Colombia/Analysis/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome \
#         -U $i \
#         -S ${iSUB}.sam ) 2>${i}.log
    
#     done

#     source activate multiqc 
#     multiqc -n bowtie2 .

# 	cd ..


# }

# alignBowtie2




dedupBAM(){

	cd primary-BAMS/dupMarked-BAMS
	# alignment stats etc. on dupMarked no MT bams
	
	for i in *.dupMarked.sorted.bam
	
	do
	
	iSUB=`echo $i | cut -d "." -f1`
	samtools view -b -h -F 0X400 $i > ${iSUB}.DEDUP.bam
	
	done

	for i in *.DEDUP.bam; do samtools index $i ; samtools idxstats $i > `echo $i | cut -d "." -f1`.DEDUP.idxstats; done
	for i in *.DEDUP.bam; do samtools flagstat $i > `echo $i | cut -d "." -f1`.DEDUP.flagstat; done

	
	mkdir dedup-BAMS
	mv *.DEDUP* dedup-BAMS/
	mv dedup-BAMS ../../
	
	cd ../..

}


# macs2 callpeak -t 2815-C.DEDUP.bam \
# 	-c 2815-I.DEDUP.bam \
#  	-f BAM -g 2.7e9 \
# 	-n 2815-C \
# 	--outdir macs2 2> macs2/2815-C-macs2.log







##### FUNC CALLS #################################################################################

#trim
#align
#sort
# markDups
 dedupBAM