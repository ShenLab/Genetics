#!/bin/bash

ANNOVAR='/ifs/data/c2b2/ys_lab/wm2313/softwares/annovar'
$ANNOVAR/table_annovar.pl $1 $ANNOVAR/humandb/ -buildver hg19 -out $2 -remove -protocol refGene,esp6500si_all,1000g2012apr_all,ljb_all -operation g,f,f,f -nastring NA -csvout -otherinfo

