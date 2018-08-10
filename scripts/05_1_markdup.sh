#!/bin/bash

data=		#path to directory containing bam files (resulting from step 04)
script=		#path to 05_2_markdup.qsub file 

for file in $(ls ${data}/*F4_F256_q5_f2_sorted.bam)
do
	base=${file##*/}
	name=${base%%_*}
        toEval="cat ${script} | sed 's/__BASE__/$base/g'"; eval $toEval > ${script%.*}_$name.qsub 
	qsub ${script%.*}_$name.qsub
done
