#!/bin/bash


bam=$1
REF=$2

if [[ $REF == "" ]]; then
    source $HOME/code/Genetics/pipeline_settings_hg19.sh
fi

if [[ ! -e $bam.flagstat ]]; then
	samtools flagstat $bam > $bam.flagstat
fi	

cram=`echo $bam | rev |cut -d '.' -f2-  | rev`".cram"

crumble -9 -O cram,reference=$REF $bam $cram
samtools index $cram
samtools flagstat $cram > $cram.flagstat

g=`diff $bam.flagstat $cram.flagstat`
echo $g
