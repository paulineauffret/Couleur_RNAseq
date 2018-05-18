#!/bin/bash

WD=			#path to filtered bam files directory (resulting from step 03)
OUT=			#path to ouptut directory
SCRIPT=                 #path to scripts directory
HEADER=                 #path to header.txt file (containing PBS parameters)
GFF=			#path to gff file
EXT=			#filtered bam extension
SAMTOOLS_ENV=   	#path to samtools version 1.4.1 conda environment if needed
SAMTOOLS=       	#path to samtools version 1.4.1
HTSEQ_ENV=		#path to htseq version 0.6.1 conda environment if needed
HTSEQ_COUNT=		#path to htseq-count version 0.6.1

cd ${WD}

mkdir -p ${OUT}

for file in $(ls *${EXT})
do
	cp ${HEADER} ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "cd ${WD}" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "${SAMTOOLS_ENV}" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "${SAMTOOLS} sort ${file} > ${file%.*}_sorted.bam" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "${SAMTOOLS} index ${file%.*}_sorted.bam > ${file%.*}_sorted.bam.bai" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "${HTSEQ}" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	echo "${HTSEQ_COUNT} -f 'bam' -s 'no' -r 'pos' -t 'CDS' -i 'Name' ${file%.*}_sorted.bam ${GFF} >> ${OUT}/${file}_htseq-count.txt" >> ${SCRIPT}/counting_htseq_${file}.qsub ;
	qsub ${SCRIPT}/counting_htseq_${file}.qsub ;
done

