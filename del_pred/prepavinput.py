## Usage: python prepavinput.py <eval dataset csv> > <eval dataset avinput>
## Purpose: reformats variant datasets into ANNOVAR input format
## Output: <eval dataset .avinput> -- chr#, start, end, ref, alt, comments
import sys

f = open(sys.argv[1],'r')

for line in f:
	tmp = line.strip().split(',')
	chr = tmp[0][3:]
	start, end = tmp[1], tmp[1]
	ref, alt = tmp[2], tmp[3]
	lab = tmp[4]
	comments = "|".join(tmp[5:])
	out = [chr, start, end, ref, alt, lab, comments]
	print "\t".join(out)
f.close()
