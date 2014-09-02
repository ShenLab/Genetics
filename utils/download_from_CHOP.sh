#!/bin/bash
#$ -cwd

list=$1
user=$2

for f in `cat $list`
do 
g="/pcgc"$f

movedat -a 3 -b 1048576 -K $user@resrhdxp01.research.chop.edu:$g ./

done

