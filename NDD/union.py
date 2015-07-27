## Usage: python union.py <list1> <list2> <list3> > <union list>
## Purpose: union operation on 3 lists, removes duplicates
## Output: union of 3 lists
import sys

f1 = open(sys.argv[1],'r')
f2 = open(sys.argv[2],'r')
f3 = open(sys.argv[3],'r')

l1 = []
for line in f1:
    l1.append(line.strip().split()[0])
f1.close()

l2 = []
for line in f2:
    l2.append(line.strip().split()[0])
f2.close()

l3 = []
for line in f3:
    l3.append(line.strip().split()[0])
f3.close()

full = l1+l2+l3
uniq = list(set(full))
for line in uniq:
    print line.strip()
