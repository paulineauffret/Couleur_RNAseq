#!/bin/bash

#Global variables
WD=		#path to working directory containing bam files (resulting from step 02)
SCRIPT=		#path to scripts directory
HEADER=         #path to header.txt file (containing PBS parameters)
SAMTOOLS_ENV=	#path to samtools version 1.4.1 conda environment if needed
SAMTOOLS=	#path to samtools version 1.4.1
SAMTOOLS_FILTERING_PARAMS="-F 4 -F 256 -q 5 -f2"	#samtools filtering parameters
EXT="F4_F256_q5_f2"					#filtered bam file extension


cd $WD

mkdir -p ${WD}/filtered_bam

for file in $(ls *.bam)
do
	cp ${HEADER} ${SCRIPT}/filter_${file}.qsub ;
	echo "${SAMTOOLS_ENV}" >> ${SCRIPT}/filter_${file}.qsub ;
	echo "${SAMTOOLS} view -b ${SAMTOOLS_FILTERING_PARAMS} ${WD}/${file} > ${WD}/${file%.*}_${EXT}.bam " >> ${SCRIPT}/filter_${file}.qsub ;
	echo "${SAMTOOLS} flagstat ${WD}/${file%.*}_${EXT}.bam > ${WD}/${file%.*}_${EXT}.flagstat " >> ${SCRIPT}/filter_${file}.qsub ;
	echo "${SAMTOOLS} flagstat ${WD}/${file} > ${WD}/${file}.flagstat " >> ${SCRIPT}/filter_${file}.qsub ;
	qsub ${SCRIPT}/filter_${file}.qsub ;
done

