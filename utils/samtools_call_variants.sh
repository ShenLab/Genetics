## call mutations from BAM file using samtools

samtools="samtools_1.3"
bcftools="bcftools"
#location="13:32,889,617-32,973,805"
location="12:24,386,753-26,403,863"
ref="genome.fa"
bamlist=$1
##"KC367_human.bam"

$samtools mpileup -r $location -Bug -t DP,ADF,ADR,SP,AD  -d 100000 -m 1 -F 0.0001  -R -f $ref $bamlist | $bcftools call -O v -v -c -p 1.1  > $bamlist.vcf 
###$bcftools view  $bamlist.bcf > $bamlist.vcf

