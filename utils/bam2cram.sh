#!/bin/bash

source $HOME/code/Genetics/pipeline_settings_hg19.sh

bam=$1
cram=`echo $bam | rev |cut -d '.' -f2-  | rev`".cram"

echo "convert $bam to $cram using samtools"
samtools view -T $REF -C -o $cram $bam

echo "make index of $cram"
samtools index $cram
 
