## Usage: python intersect.py <list1> <list2> > <overlap list>
## Purpose: intersects two lists
## Output: prints intersection of two lists

import sys

f1 = open(sys.argv[1],'r')
f2 = open(sys.argv[2],'r')

l1 = []
for line in f1:
    l1.append(line.strip())
f1.close()

l2 = []
for line in f2:
    l2.append(line.strip().split()[0])
f2.close()

inter1 = list(set(l2).intersection(l1))
for line in inter1:
    print line.strip()
