#!/bin/bash
#$ -cwd

dir=$1 ## the dir of fastq files
## /ifs/home/c2b2/ys_lab/wm2313/gtex/Heart_reads/fastq/
for r1 in $dir/*_1.fastq.gz; do
    r1_base=`basename $r1 _1.fastq.gz`
    r2=$dir/${r1_base}_2.fastq.gz
    
    if [[ -e $r2 ]]; then
       echo $r1
       echo $r2
	   sh handleStar.sh $r1 $r2;
    fi
done
