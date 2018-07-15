#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=5000::
#$ -l h_vmem=3G

list=$1

for f in `cat $list`
do
if [[ -e $f ]]; then
bash bam2cram_hg38.sh $f 	
echo "$f is done"

echo "" > $f
echo "" > $f.bai 
fi 
done
