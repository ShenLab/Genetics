#!/bin/bash
#$ -cwd

USAGE="Usage: $0 -i <input BAM> -o <output BAM> \n"


while getopts i:o:h opt
do case "$opt" in
	i)      input="$OPTARG";;
	o)      output="$OPTARG";;
	h)      echo $USAGE
	    exit 1;;
    esac
done

export PATH=$PATH:$HOME/data/usr/bin/

if [[ $output == "" ]]; then
    output="markDup.bam"
fi

if [[ $input != "" ]]; then
    samtools view -b $input | bammarkduplicates2    O=$output
else 
    bammarkduplicates2 I=$input O=$output
fi

samtools index $output
samtools flagstat $output > $output.flagstat


