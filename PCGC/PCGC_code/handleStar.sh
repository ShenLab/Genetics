#!/bin/bash -l
#$-cwd

r1=$1
r2=$2
code=/ifs/home/c2b2/ys_lab/wm2313/code/
log=starlogs

if [ ! -d  $log ] ; then mkdir -p $log; fi

r1_base=`basename $r1 .gz`
r2_base=`basename $r2 .gz`

echo $r1_base

## qsub -l mem=11G,time=16:: -P shen-chung -pe smp 4 -m e -o $log/star.$r1_base.o -e $log/star.$r1_base.e  -N star.$r1_base $code/star.sh $r1 $r2
qsub -l mem=6G,time=14:: -P ys_lab -pe smp 6 -o $log/star.$r1_base.o -e $log/star.$r1_base.e  -N star.$r1_base $code/star.sh $r1 $r2
