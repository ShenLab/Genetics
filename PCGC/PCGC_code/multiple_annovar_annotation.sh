#!/bin/bash

FILES='/ifs/scratch/c2b2/ys_lab/hq2130/beforannovar/*.avinput'
OUT='/ifs/scratch/c2b2/ys_lab/hq2130/afterannovar'

for f in $FILES
do

af=$OUT${f:44:15}

echo $af
echo $f

qsub  -l mem=8G,time=4:: -S /bin/bash -cwd -N my_annovar -pe orte 4 ./annovar_annotate.sh $f $af


done
