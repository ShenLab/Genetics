## Usage: python union.py <list1> <list2> <fdr cutoff> > <union list>
## Purpose: union operation on 2 gene lists, removes duplicates and entries with fdr above cutoff value
## Output: union of 2 lists
import sys

f1 = open(sys.argv[1],'r')
f2 = open(sys.argv[2],'r')
fdrcut = float(sys.argv[3])

f1l = []
for line in f1:
    tmp = line.strip().split()
    f1l.append(tmp[0])
f1.close()

f2l = []
f2d = dict()
for line in f2:
    tmp = line.strip().split()
    if not tmp[0] == "gene":
        f2l.append(tmp[0])
        fdr = float(tmp[12])
        f2d[tmp[0]] = fdr
f2.close()

inter = list(set(f1l+f2l))
yesfdr = 0
for g in inter:
    if g in f2d:
        if f2d[g] <= fdrcut:
            print g+","+str(f2d[g])
