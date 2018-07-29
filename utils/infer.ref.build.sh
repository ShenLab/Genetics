#!/bin/bash
#$ -cwd

list=$1

for f in `cat $list`
do 

g=`samtools view -h $f | head -200 | egrep "^\@SQ" | egrep -w "LN\:248956422"`
h=`samtools view -h $f | head -200 | egrep "^\@SQ" | egrep -w "LN\:249250621"`
 if [[ $g =~ \@SQ*  ]]; then
     ref="hg38"
elif [[ $h =~ \@SQ* ]]; then
     ref="hg19"
else
     ref="bad_file"
fi

 echo "$f $ref"

done

