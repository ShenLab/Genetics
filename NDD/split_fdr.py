## Usage: python split_fdr.py external_DD.txt <lower bound> <upper bound>
## Purpose: splits external_DD.txt genes file into a smaller gene list based on FDR cutoffs
## Output: prints list of gene symbols with FDR falling in the specified range

import sys

f = open(sys.argv[1],'r')

fdrs = []
for line in f:
    tmp = line.strip().split()
    if not tmp[0] == "gene":
        fdrs.append(float(tmp[12]))
        if float(sys.argv[2]) < float(tmp[12]) <= float(sys.argv[3]):
            print tmp[0]
f.close()
