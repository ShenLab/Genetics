#!/bin/bash

bam=$1

samtools view -H $bam | egrep -w "^\@SQ\s+SN" > $bam.refinfo 

snline=`wc -l $bam.refinfo | awk '{print $1}'`
hs37d5=`grep -c hs37d5 $bam.refinfo`
chr=`grep -c chr $bam.refinfo`
ebv=`grep -c NC_007605 $bam.refinfo`
chr1hg19=`grep -c LN\:249250621 $bam.refinfo`
# echo $chr1hg19
# echo $chr

if [[ $chr1hg19 == "1" && $chr == "0" ]]; then
    echo "$bam    hs37d5"   ### hs37d5, decoy patch for hg19, used by RGN
elif [[ $chr1hg19 == "1" && $chr != "0" ]]; then
    echo "$bam    ucsc_hg19"
else
    echo "$bam    other"
fi

