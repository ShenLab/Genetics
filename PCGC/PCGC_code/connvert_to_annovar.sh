#!/bin/bash
Files="/Users/hongjian/Documents/research/data/*.vcf"
Out="/Users/hongjian/Documents/research/annovar/"
for f in $Files

do
 af=$Out${f:40:14}'.avinput'
 echo $af

 convert2annovar.pl -format vcf4old $f -outfile $af -include

done