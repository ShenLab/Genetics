#!/bin/bash -l
#$ -cwd

bam=$1
gtf=$2

gsnap_dir=/ifs/home/c2b2/ys_lab/wm2313/data/GSNAP
reference=$gsnap_dir/hg19.fa
bam_base=`basename $bam .bam`
sample=${bam_base}
reheaderbam=${bam_base}".reheader.bam"

if [[ ! -e $reheaderbam.bai ]]; then

echo reheader the bam files for GATK...
~/data/softwares/jdk1.7.0_60/bin/java -jar -Xmx3g /ifs/home/c2b2/ys_lab/wm2313/data/softwares/picard/AddOrReplaceReadGroups.jar I=$bam O=$reheaderbam LB=XX PL=XX PU=XX SM=XX
echo index the new bam file..
~/data/softwares/samtools-0.1.19/samtools index $reheaderbam

fi


echo perform RNA-SeQC..
~/data/softwares/jdk1.7.0_60/bin/java -jar -Xmx10g /ifs/home/c2b2/ys_lab/wm2313/data/softwares/RNA-SeQC_v1.1.7.jar -o ./"RNA-SeQC" -gatkFlags "-DBQ 1" -r $reference -t $gtf  -s "$sample|$reheaderbam|NA"



## get a summary of the mapping quality
cd RNA-SeQC

summary=../../rnaseqc.info
if [[ ! -e $summary ]]; then
     echo -e "sample\tReadsNO\tmappingRate\tintraRate\texonicRate\tintronicRate\tintergenicRate\treadsLength" > $summary
fi

ID=`echo ${bam_base} | cut -d '.' -f1-3`

if [[ -e metrics.tsv ]]; then

    data=`tail -n +2 metrics.tsv`
    intraRate=`echo $data | cut -d " " -f5,5`
    exonicRate=`echo $data | cut -d " " -f7,7`

    mappingRate=`echo $data | cut -d " " -f8,8`

    readsLength=`echo $data | cut -d " " -f12,12`
    mappedNO=`echo $data | cut -d " " -f17,17`

    intergenicRate=`echo $data | cut -d " " -f18,18`
    intronicRate=`echo $data | cut -d " " -f28,28`
    readsNO=`echo "$mappedNO/$mappingRate" | bc`

    echo -e "$ID\t$readsNO\t$mappingRate\t$intraRate\t$exonicRate\t$intronicRate\t$intergenicRate\t$readsLength" >> $summary
    
    cd ..
    rm -f $reheaderbam
    rm -f $reheaderbam.bai

else
    echo -e "$ID\tfails on RNA-SeQC" >> $summary
fi
