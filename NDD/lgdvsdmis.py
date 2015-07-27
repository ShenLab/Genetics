## Usage: python predict5.py cases_metaSVM.csv phenotypecombined.csv fdrbins/fdr_to0.5.txt
## Purpose: Look at prevalence of different mutation types in risk genes list for patients with both certain and uncertain diagnosis
## Output: Among certain/uncertain diagnosis patients, prints counts for different mutation types within each group
import sys

f1 = open(sys.argv[1],'r') # mutations file
f2 = open(sys.argv[2],'r') # clinical phenotypes file
f3 = open(sys.argv[3],'r') # risk genes list

# read in risk genes
tmp = []
for line in f3:
    tmp.append(line.strip().split()[0])
f3.close()
genes = list(set(tmp)) # delete duplicates

# read in phenotypes file
phend = dict()
for line in f2:
    tmp = line.strip().split(",")
    if not tmp[0] == "ID":
        if tmp[3] == "NDD":
            ndd = 1
        elif tmp[3] == "non-NDD":
            ndd = 0
        elif tmp[3] == "uncertain":
            ndd = -1
        if tmp[4] == "Extra-pos":
            ec = 1
        elif tmp[4] == "Extra-neg":
            ec = 0
        phend[tmp[0]] = [ec, ndd]
f2.close()

# read in mutations file
idx = dict()
mutd = dict()
anno = dict()
for line in f1:
    tmp = line.strip().split(',')
    if tmp[0] == "Subject": # header line
        for i in xrange(len(tmp)):
            idx[tmp[i]] = i
    else: # data lines
        if not tmp[0] in mutd:
            mutd[tmp[0]] = [[],[],[],[]] # LGD, Dmis, NDmis, other
        if not tmp[0] in anno:
            anno[tmp[0]] = []
        ec = phend[tmp[0]][0]
        ndd = phend[tmp[0]][1]
        mutgene = tmp[5]
        type = ''
        if tmp[idx["ExonicFunc.refgene"]] in ['stopgain','stoploss','frameshift insertion','frameshift deletion'] or 'splicing' in tmp[idx["Func.refgene"]]:
            type = 'LoF'
            mutd[tmp[0]][0].append(mutgene)
        elif tmp[idx["ExonicFunc.refgene"]] in ['nonsynonymous SNV']:
            if tmp[idx["MetaSVM_pred"]] == "D": # Dmis
                type = 'Dmis'
                mutd[tmp[0]][1].append(mutgene)
            else:
                type = 'NDmis'
                mutd[tmp[0]][2].append(mutgene)
        else:
            type = 'other'
            mutd[tmp[0]][3].append(mutgene)
        anno[tmp[0]] = [mutd[tmp[0]], ec, ndd]
f1.close()

lofct, dmisct, ndmisct, otherct = 0, 0, 0, 0
certaincts = [0,0,0,0]
uncertaincts = [0,0,0,0]
certains = 0
uncertains = 0
for pt in anno:
    nddstat = anno[pt][2]
    if nddstat == -1:
        uncertains += 1
    elif nddstat in [0,1]:
        certains += 1
    ecstat = anno[pt][1]
    lofgenes, dmisgenes, ndmisgenes, othergenes = anno[pt][0][0], anno[pt][0][1], anno[pt][0][2], anno[pt][0][3]
    if len(lofgenes) > 0:
        if len(list(set(lofgenes).intersection(genes))) > 0:
            lofct += 1
            if nddstat == -1:
                uncertaincts[0] += 1
            elif nddstat in [0, 1]:
                certaincts[0] += 1
    if len(dmisgenes) > 0:
        if len(list(set(dmisgenes).intersection(genes))) > 0:
            dmisct += 1
            if nddstat == -1:
                uncertaincts[1] += 1
            elif nddstat in [0, 1]:
                certaincts[1] += 1
    if len(ndmisgenes) > 0:
        if len(list(set(ndmisgenes).intersection(genes))) > 0:
            ndmisct += 1
            if nddstat == -1:
                uncertaincts[2] += 1
            elif nddstat in [0, 1]:
                certaincts[2] += 1
    if len(othergenes) > 0:
        if len(list(set(othergenes).intersection(genes))) > 0:
            otherct += 1
            if nddstat == -1:
                uncertaincts[3] += 1
            elif nddstat in [0, 1]:
                certaincts[3] += 1
print "certain cases: "+str(certains)
print "uncertain cases: "+str(uncertains)
print "total:"
print [lofct, dmisct, ndmisct, otherct]
print "certain cases:"
print certaincts
print "uncertain cases:"
print uncertaincts


