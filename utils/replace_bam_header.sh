#!/bin/bash
#$ -cwd

export PATH=$PATH:$HOME/usr/bin/

header=$1
bam=$2

samtools reheader $header $bam > $bam.temp
mv $bam.temp $bam
samtools index $bam 

