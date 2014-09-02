#!/bin/bash
#$ -cwd

USAGE="Usage: $0 -i <input BAM> -o <output BAM> [-t n_threads]\n"

threads=1

while getopts i:o:t:h opt
do case "$opt" in
	i)      input="$OPTARG";;
	o)      output="$OPTARG";;
	t)      threads="$OPTARG";;
	h)      echo $USAGE
	    exit 1;;
    esac
done

export PATH=$PATH:$HOME/data/usr/bin/

if [[ $output == "" ]]; then
    output="markDup.bam"
fi

if [[ $input != "" ]]; then
    samtools view -b $input | bammarkduplicates2    O=$output markthreads=$threads
else 
    bammarkduplicates2 I=$input O=$output markthreads=$threads
fi

samtools index $output
samtools flagstat $output > $output.flagstat


