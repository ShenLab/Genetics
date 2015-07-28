## Usage: python union.py list1.txt list2.txt .... > unionlist.txt
## Purpose: union operation on a variable number of list files, removes duplicates
## Output: union of all argument lists
import sys

tmp = []

for arg in sys.argv[1:]:
	f = open(arg,'r')
	for line in f:
		tmp1 = line.strip().split()
		tmp.append(tmp1[0])
	f.close()

uniq = list(set(tmp))
for l in uniq:
	print l
