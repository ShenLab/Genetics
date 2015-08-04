## Usage: python getnon.py metasvm_trainingdataset.genes.csv eval_domains_training.csv > non_domains_training.csv
## Purpose: Generates set of variants that do NOT have annotated domain/3D structure information
## Output: chr,pos,ref,alt,label,gene

import sys

metaf = open(sys.argv[1],'r')
annof = open(sys.argv[2],'r')

present = []
for line in annof:
	tmp = line.strip().split(',')
	key = tmp[0]+tmp[1]+tmp[2]+tmp[3] #chr,pos,ref,alt
	present.append(key)
annof.close()

for line in metaf:
	tmp = line.strip().split('\t')
	check = tmp[0]+tmp[1]+tmp[2]+tmp[3]
	if not check in present:
		print ','.join(tmp)
metaf.close()