#!/bin/bash
#$ -cwd

export PATH=$PATH:$HOME/usr/bin/

bam=$1

setting=$2

. $setting

bedtools multicov -bams $bam -bed $ExonFile -q 20 > $bam.canoes.reads.txt
