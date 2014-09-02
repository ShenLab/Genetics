#!/bin/bash
#$ -cwd

export PATH=$PATH:$HOME/usr/bin/

fq=$1

seqtk seq -Q64 -V $fq  | gzip > $fq.new.gz

mv $fq.new.gz $fq

