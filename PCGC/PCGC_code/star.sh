#!/bin/bash -l
#$ -cwd

R1=$1
R2=$2

## check the length of reads to choose star configuration
read=`zcat -d $R1| head | tail -n +10`
readLength=`echo ${#read}`
echo "fastq read length: "$readLength

out="StarBam"
if [ ! -d $out ] ; then mkdir -p $out; fi

R1_base=`basename $R1`
cd $out

if [ -d ${R1_base} ] ; then rm -r ${R1_base}; fi

mkdir -p ${R1_base}
cd ${R1_base}

echo "star align begin..."
date

if [[ $readLength -gt 60 ]];
    then
    genomeDB="~/scratch/softwareMakeup/genome100"
    seedStart=50
    else
    genomeDB="/ifs/home/c2b2/ys_lab/wm2313/data/softwares/STAR_2.3.0e/genome"
    seedStart=30
fi

/ifs/home/c2b2/ys_lab/wm2313/data/softwares/STAR_2.3.0e/STAR --genomeDir $genomeDB --sjdbGTFfile /ifs/home/c2b2/ys_lab/wm2313/data/GSNAP/chr_genes.gtf --readFilesCommand zcat --seedSearchStartLmax $seedStart --runThreadN 4 --outSAMunmapped Within  --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical  --readFilesIn $R1 $R2 --outStd SAM --outSAMmode Full | ~/data/softwares/samtools-0.1.19/samtools view -bS -> ${R1_base}.bam

## --outSAMstrandField intronMotif is used to add XS strand attribute for alignments that contain splice junctions.
## --outFilterIntronMotifs RemoveNoncanonical filter out alignments that contain non-canonical junctions

echo "star alignment finished: "
date

bam=${R1_base}.bam
sortbam=${R1_base}.hg19UCSC.star.sort.bam
~/data/softwares/samtools-0.1.19/samtools sort $bam ${R1_base}.hg19UCSC.star.sort
echo "samtools sort finished: "
date

if [[ -e $sortbam ]]; 
then
rm -f $bam
fi

~/data/softwares/samtools-0.1.19/samtools index $sortbam
echo "samtools index finished: "
date

if [[ -e ${R1_base}".hg19UCSC.star.sort.bam.bai" ]];

    then
        echo ${R1_base} done! >> ../starbam.checklist

    else
        echo echo ${R1_base} failed! >> ../starbam.checklist
fi



logs="analysis_logs"
if [ ! -d  $logs ] ; then mkdir -p $logs; fi
currentDir=`pwd`
parentDir=`dirname $currentDir`
## gtf=/ifs/home/c2b2/ys_lab/wm2313/data/GSNAP/chr_genes.gtf
gtf=/ifs/home/c2b2/ys_lab/wm2313/data/GSNAP/gencode.v19.annotation.gtf

## peform alignment QC
qsub -l mem=5G,time=5:: -pe smp 6 -o $logs/${R1_base}.RNAseqc.o -e $logs/${R1_base}.RNAseqc.e -N RNASEQC.${R1_base} ~/code/star_RNA-SeQC_pipe.sh $sortbam $gtf

## compute FPKM with cufflinks
mkdir $parentDir/fpkms ## where to put the fpkm files(gene level)
qsub -l mem=1G,time=16:: -pe smp 6 -m e -o $logs/${R1_base}.cufflink.o -e $logs/${R1_base}.cufflink.e -N cufflinks.${R1_base} ~/code/cufflinks.sh $sortbam $gtf $parentDir/fpkms

