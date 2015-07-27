## Usage: python parsesex.py gender_phenotypes.ped phenotypecombined.csv
## Purpose: Prints summary stats for number of males and females, and for each, the number of NDD/nonNDD/uncertain
## Output: prints summary stats
import sys

f1 = open(sys.argv[1],'r')
f2 = open(sys.argv[2],'r')

gend = dict()
for line in f1:
    tmp = line.strip().split()
    if len(tmp[2]) > 2:
        if tmp[4] == "1":
            s = "M"
        elif tmp[4] == "2":
            s = "F"
        if not tmp[1] in gend:
            gend[tmp[1]] = s
f1.close()

gendct = 0
uncct = 0
ndd = [[],[]]
non = [[],[]]
unc = [[],[]]
for line in f2:
    tmp = line.strip().split(',')
    if not tmp[0] == "ID":
        if tmp[0] in gend:
            gendct += 1
            if tmp[3] == "NDD":
                if gend[tmp[0]] == "M":
                    ndd[0].append(tmp[0])
                elif gend[tmp[0]] == "F":
                    ndd[1].append(tmp[0])
            elif tmp[3] == "non-NDD":
                if gend[tmp[0]] == "M":
                    non[0].append(tmp[0])
                elif gend[tmp[0]] == "F":
                    non[1].append(tmp[0])
            elif tmp[3] == "uncertain":
                if gend[tmp[0]] == "M":
                    unc[0].append(tmp[0])
                elif gend[tmp[0]] == "F":
                    unc[1].append(tmp[0])
f2.close()
print "total: "+str(gendct)
print "male:female"
print str(len(ndd[0])) + " : " + str(len(ndd[1]))
print str(len(non[0])) + " : " + str(len(non[1]))
print str(len(unc[0])) + " : " + str(len(unc[1]))
