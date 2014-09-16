#!/bin/bash
#$ -cwd

httplink=$1
wget -N -r -nH --cut-dirs=2 -np -l 2 -A gz --execute robots=off $httplink
