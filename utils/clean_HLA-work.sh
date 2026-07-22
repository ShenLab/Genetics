### clean working dir for HLA allele inference by WGS data

## files: 
# to rm:  extraction*     *.fastq            

## files to keep: remapped_with_a.bam.bai    hla  
##             remapped_with_a.bam

# compress: reads_per_level.txt 

dir=$1

if [[ -e $dir/remapped_with_a.bam.bai && -e $dir/reads_per_level.txt  ]]; then
## clean
	rm -rf $dir/extraction* $dir/*.fastq
#	gzip $dir/reads_per_level.txt 
	echo "$dir is cleaned"
else
	echo "$dir not complete, keep intact" 
fi


