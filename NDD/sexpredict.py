## Usage: python sexpredict.py <phenotypes_wkc file> <sex/gender file>
## Purpose: for each major category of extracardiac abnormality, print fraction of male/female
## Output: for each EC major category, print number of males and females

import sys

ecf = open(sys.argv[1],'r')
sexf = open(sys.argv[2],'r')

patd = dict() # keys = patients
catd = dict() # keys = categories
for line in ecf:
    tmp = line.strip().split(',')
    if not tmp[0] == "ID":
        if not tmp[0] in patd:
            patd[tmp[0]] = []
        maj_cats = tmp[1].split(';')
        tmp1 = []
        for i in maj_cats:
            tmp2 = []
            for j in i.split('.'):
                if not j in ["General","Other","Abnormality","Abnormalities","Present"]:
                    tmp2.append(j)
                cat = '.'.join(tmp2)
                if not cat in catd:
                    catd[cat] = []
                catd[cat].append(tmp[0])
                tmp1.append(cat)
        patd[tmp[0]] = list(set(tmp1))
ecf.close()

sexd = dict()
for line in sexf:
    tmp = line.strip().split()
    sexd[tmp[0]] = tmp[1]
sexf.close()

for c in catd:
    m = 0
    f = 0
    tmp = list(set(catd[c]))
    for pat in tmp:
        if sexd[pat] == "M":
            m += 1
        elif sexd[pat] == "F":
            f += 1
    print c + " : " + str(m) + " males " + str(f) + " females "
    #print c + " : " + str(float(m)/float(m+f))+ " : "+str(float(f)/float(m+f))

