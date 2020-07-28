#!/bin/bash

list=$1
ref=$2


for f in `cat $list`
do 
	g=`bash crumble_convert.sh $f $ref`

if [[ $g == "" ]]; then
        echo "" > $bam
else
        echo "$cram is different to $bam"
fi    

done
