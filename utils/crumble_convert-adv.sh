#!/bin/bash

bam=$1
REF=$2

if [[ $REF == "" ]]; then
    source $HOME/code/Genetics/pipeline_settings_hg19.sh
fi

if [[ ! -e $bam.flagstat ]]; then
	samtools flagstat $bam > $bam.flagstat
fi	

if [[ ! -e $bam.bai ]]; then
	samtools index $bam
fi

cram=`echo $bam | rev |cut -d '.' -f2-  | rev`".cram"


if [[ -e $cram ]]; then ## cram already exist
	echo "$cram already exists" 
	exit
fi

crumble -9 -O cram,reference=$REF $bam $cram
samtools index $cram
samtools flagstat $cram > $cram.flagstat

g=`diff $bam.flagstat $cram.flagstat`
fs=`ls -l $cram.flagstat | awk '{print $5}'`
if [[ $g != "" || $fs == "0" ]]; then
    echo "$bam and $cram do not match"
else	
    echo "successfully converted bam to cram $bam"
fi
