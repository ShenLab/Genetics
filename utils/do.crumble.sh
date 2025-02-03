#!/bin/bash

list=$1
ref=$2

for f in `cat $list`
do 
    echo "working on $f"
    bash $HOME/code/Genetics/utils/crumble_convert.sh $f $ref

#    if [[ $g == "" ]]; then
#	echo "$f is converted to cram format by crumble"
#        echo "" > $f   # remove the content of the original BAM file to save disk space
#    else
#        echo "$f is kept intact because the new cram file does not have the same flagstat"
#    fi    

done
