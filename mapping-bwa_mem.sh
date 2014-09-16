#!/bin/bash
#$ -cwd
uname -a
echo "Start $0:`date`"

# default values
fastq1=""
fastq2=""

minScore=25
minSeed=19
band=150
platform="Illumina"
threads=2
sampleName=""
readgroup=""
output=""
chain="0"  # whether call downstream analysis (default 0 means not). 

## place
groups=`groups $USER | awk '{print $3}' | tr '[a-z]' '[A-Z]'`

### Note: bwa sample uses about 3.5G RAM

HTSCMD="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/htscmd"

export PATH=$PATH:$HOME/usr/bin/

USAGE="Usage: $0 -i foo_1.fastq  -s global_setting [ -p foo_2.fastq ] [ -k minSeedLength] [ -w  bandWidth ] [ -y ID] [ -z readgroup] [ -n sampleName] [ -f platform] [-o output.bam] \n"

while getopts i:b:p:T:k:w:t:s:y:z:n:f:o:c:A:h opt
do      
    case "$opt" in
	i) fastq1="$OPTARG";;
	b) bam="$OPTARG";;
	p) fastq2="$OPTARG";;
	T) minScore="$OPTARG";;
	k) minSeed="$OPTARG";;
	w) band="$OPTARG";;
	t) threads="$OPTARG";;
	s) setting="$OPTARG";;  
	y) ID="$OPTARG";;
	z) readgroup="$OPTARG";;
	n) sampleName="$OPTARG";;
	f) platform="$OPTARG";;
	o) output="$OPTARG";;
	c) chain="$OPTARG";;
	A) AUTO="$OPTARG";;
	h)    echo $USAGE
          exit 1;;
  esac
done


#echo $fastq1 
#echo $setting

if [[ ( $fastq1 == "" && $bam == "" ) || $setting == "" ]]; then
    echo $USAGE
    exit 1
fi

## load global setting
. $setting

bwaVersion=`$BWA 2>&1 | grep Version | awk '{print $2}'`
echo "bwa version: $bwaVersion"


if [[ $readgroup == "" ]]; then
    readgroup=`basename $fastq1  | sed 's/.gz$//'  | sed 's/.fastq$//'  | sed 's/.txt$//' | sed s'/.R1$//' | sed 's/\_R1$//'`
fi

if [[ $sampleName == "" ]]; then
    sampleName=$readgroup
fi

if [[ $ID == "" ]]; then
    ID=$readgroup
fi

if [[ $output == "" ]]; then
    output=$fastq1.sorted
fi


## read group specification:
##          -r STR   read group header line such as `@RG\tID:foo\tSM:bar' [null]
rgheader="@RG\tID:$ID\tSM:$sampleName\tLB:$readgroup\tPL:$platform\tCN:COLUMBIA_"$groups

# stream fastq
# samtools bam2fq reads.bam | bwa mem -p ref.fa -

######### bwa mem
bwacmd="$BWA mem -T $minScore -k $minSeed -w $band -t  $threads -R $rgheader -M   $REF"

if [[ $fastq1 != "" ]]; then
    cmd="$bwacmd $fastq1 $fastq2" 
    echo $cmd "|" "$SAMTOOLS view -bS - " "> $output.bam"
    $cmd  | $SAMTOOLS view -bS -   > $output.bam
elif [[ $bam != "" ]]; then
   #  cmd1="$SAMTOOLS bam2fq $bam" 
    cmd1="$HTSCMD bamshuf -uO  $bam $output.tmp" 
    cmd2="$HTSCMD bam2fq -as $output.se.fastq.gz  - "
    cmd3="$bwacmd -p - "
    echo $cmd1 "|" $cmd2 "|" $cmd3 "|" "$SAMTOOLS view -bS - " "> $output.bam"
    $cmd1 | $cmd2 | $cmd3 | $SAMTOOLS view -bS -   > $output.bam
fi


date
echo "bwa alignment complete. Sorting the bam file ..."

echo -e "@HD\tVN:1.0\tGO:none\tSO:coordinate" >  $output.bam.header

$SAMTOOLS view -H $output.bam | egrep -v '^\@HD' >> $output.bam.header
echo -e "@PG\tID:bwa\tVN:$bwaVersion\tCL:$BWA mem -T $minScore -k $minSeed -w $band -t $threads -M $REF " >> $output.bam.header


$SAMTOOLS sort -@ $threads -l 6  $output.bam $output.temp

$SAMTOOLS reheader $output.bam.header $output.temp.bam | bammarkduplicates2    O=$output.bam rewritebam=1 

# $SAMTOOLS reheader $output.bam.header $output.temp.bam > $output.bam 
 
# rm -f $output.temp

# $SAMTOOLS reheader $output.header $output.sort.bam > $output.bam

# mv $output.temp $output
rm -f $output.temp.bam #  $output.sort.bam 

# index
$SAMTOOLS index $output.bam

$SAMTOOLS flagstat $output.bam > $output.bam.flagstat
$SAMTOOLS idxstats $output.bam > $output.bam.idxstats
echo "$0 done"
date

