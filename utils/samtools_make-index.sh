#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=2::
#$ -l h_vmem=1G

cram=$1

samtools index $cram 

