#!/bin/bash
#$ -cwd
bam=$1
gtf=$2
dir=$3
out="Cufflinksout"

if [ ! -d $out ]; then mkdir -p $out; fi
if [ ! -d $dir ]; then mkdir -p $dir; fi

bam_base=`basename $bam .bam`
individualOut=$bam_base"_cufflinks"
cmd="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/cufflinks -o $out/$individualOut --compatible-hits-norm --GTF  $gtf $bam"
echo -e "do cufflinks with ref genes: \n $cmd"
$cmd
cp ./$out/$individualOut/genes.fpkm_tracking $dir/${bam_base}.genes.fpkm_tracking

