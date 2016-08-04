#!/bin/bash
#$ -cwd

list=$1

userID=ys2411
port="33001"

for f in `cat $list`
do 
	g=`basename $f`
	h=`echo $f | sed 's/\.\///'`
	if [[ -e $g ]]; then
		echo "$g is done"
	else
		ascp -i ~/.ssh/id_rsa -p -v -l 400m -P $port -d $userID@hgsc-aspera1.hgsc.bcm.edu:gmkf/seidman/outbound/$h . 
		echo "$g is done"
	fi
done

