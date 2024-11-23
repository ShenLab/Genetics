## purpose: clean up SPARK exome data downloaded from DNAnexus.
## 
##  most of the exome data downloaded from DNAnexus are included in iWES releases
##  this script look up the corresponding cram files in iWES and DNAnexus releases 
##  and remove the ones in DNAnexus release fold and place a softlink in place



list=$1   ## list of SPID

iWESCrams=$2
nexusCrams=$3

for f in `cat $list`
do 
 f1=`grep $f $iWESCrams`
 f2=`grep $f $nexusCrams`

 if [[ $f1 != "" && $f2 != "" ]]; then
    echo $f1 $f2

    ## replace by a link
   #  rm -rf $f2
   #  ln -s $f1 $f2 
   ### g=`echo $f2 | sed 's/cram$/crai/'`
   #  rm -rf $g
   #  ln -s $f1.crai $f2.crai
    
 fi	

done

