## Usage: python ranksplit.py external_DD.txt <lower bound> <upper bound>
## Purpose: splits external_DD.txt gene list into smaller lists based on p.min rank
## Output: genes that fall within rank range

import sys

f = open(sys.argv[1],'r')

genes = []
pmins = []
for line in f:
	tmp = line.strip().split()
	if not tmp[0] == "gene":
		pmin = float(tmp[11]) # convert to float
		tmp[11] = pmin
		genes.append(tmp)
		pmins.append(tmp[11])
f.close()
genes_sorted = sorted(genes, key=lambda x: x[11])

i = 0
for gene in genes_sorted:
	i += 1
	if float(sys.argv[2]) <= i <= float(sys.argv[3]):
		print gene[0]