## Usage: python processextracardiac.py phenotypes_wkc.csv phenotypecombined.csv
## Purpose: Process the list of extracardiac manifestations
## Output: For each EC major category, output counts for NDD/non-NDD/uncertain and posterior probability of NDD (beta binomial approximation)

import sys
from collections import Counter

f1 = open(sys.argv[1],'r') # new phenotypes file
f2 = open(sys.argv[2],'r') # old phenotypes master list
f3 = open(sys.argv[3],'r') # genders file

sex_d = dict()
for line in f3:
    tmp = line.strip().split()
    if not tmp[0] in sex_d:
        sex_d[tmp[0]] = tmp[1]
f3.close()


ndd_d = dict()
for line in f2:
    tmp = line.strip().split(',')
    ndd = 0
    if not tmp[0] == "ID":
        if tmp[3] == "NDD":
            ndd = 1
        elif tmp[3] == "non-NDD":
            ndd = 0
        elif tmp[3] == "uncertain":
            ndd = -1
        ndd_d[tmp[0]] = ndd
f2.close()

majcat = dict()
patd = dict()
for line in f1:
    tmp = line.strip().split(',')
    if not tmp[0] == "ID":
        nddstat = ndd_d[tmp[0]]
        if not tmp[0] in patd:
            patd[tmp[0]] = [[],[],nddstat]
        maj_cats = tmp[1].split(';')
        tmp1 = []
        for i in maj_cats:
            tmp2 = []
            for j in i.split('.'):
                if not j in ["General","Other", "Abnormality", "Abnormalities", "Present"]:
                    tmp2.append(j)
            cat = ".".join(tmp2)
            if not cat in majcat:
                majcat[cat] = [[],[],[]]
            if nddstat == 1: majcat[cat][0].append(tmp[0])
            elif nddstat == 0: majcat[cat][1].append(tmp[0])
            elif nddstat == -1: majcat[cat][2].append(tmp[0])
            tmp1.append(cat)
        patd[tmp[0]][0] = list(set(tmp1))
        patd[tmp[0]][1] = list(set(tmp2))
f1.close()

head = ["#EC.MajCat","NDD", "non-NDD", "uncertain", "posterior", "male_ndd", "female_ndd", "male_non", "female_non", "male_unc", "female_unc"]
print ','.join(head)
outlist = []
for i in majcat:
    # Beta-Binomial distribution parameters for estimating p
    a = 1
    b = 1
    nddl = list(set(majcat[i][0]))
    m_nddl = [pt for pt in nddl if sex_d[pt] == "M"]
    f_nddl = [pt for pt in nddl if sex_d[pt] == "F"]
    nonl = list(set(majcat[i][1]))
    m_nonl = [pt for pt in nonl if sex_d[pt] == "M"]
    f_nonl = [pt for pt in nonl if sex_d[pt] == "F"]
    uncl = list(set(majcat[i][2]))
    m_uncl = [pt for pt in uncl if sex_d[pt] == "M"]
    f_uncl = [pt for pt in uncl if sex_d[pt] == "F"]
    outlist.append([i, str(len(nddl)), str(len(nonl)), str(len(uncl)), str(float(len(nddl)+a)/float(len(nddl)+len(nonl)+a+b)), str(len(m_nddl)), str(len(f_nddl)), str(len(m_nonl)), str(len(f_nonl)), str(len(m_uncl)), str(len(f_uncl))])
    #print i+","+str(majcat[i][0])+","+str(majcat[i][1])+","+str(majcat[i][2])+","+str(float(majcat[i][0])/float(majcat[i][0]+majcat[i][1]))
for k in sorted(outlist, key=lambda x: x[4], reverse=True):
    print ','.join(k)
