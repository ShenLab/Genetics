## Usage: python predict5.py cases_metaSVM.csv phenotypecombined.csv phenotypes_extracardiac_mb_wkc.csv gender_phenotypes.parsed.txt <Gene list> 
## Purpose: for each subject, counts LoF and Mis3 mutations in risk genes list
##          then evaluates mutation counts as predictors of NDD status
## Output: counts for number of NDD, non-NDD, uncertain diagnosis patients; also outputs sex information
import sys

# read in risk genes
genef = open(sys.argv[5],'r')
tmp = []
for line in genef:
    if not line.startswith("#"):
        tmp.append(line.strip().split()[0])
genef.close()
genes = list(set(tmp)) # delete duplicate entries

# read in phenotypes file
phenf = open(sys.argv[2],'r')
phend = dict()
for line in phenf:
    tmp = line.strip().split(",")
    if not tmp[0] == "ID":
        # capture NDD diagnosis
        if tmp[3] == "NDD": 
            ndd = 1
        elif tmp[3] == "non-NDD": 
            ndd = 0
        elif tmp[3] == "uncertain": 
            ndd = -1
        # capture extracardiac status
        if tmp[4] == "Extra-pos":
            ec = 1
        elif tmp[4] == "Extra-neg": 
            ec = 0
        phend[tmp[0]] = [ec, ndd]
phenf.close()

# read in phenotypes major category file
catf = open(sys.argv[3],'r')
catd = dict()
for line in catf:
    tmp = line.strip().split(',')
    if not tmp[0] == "ID":
        nddstat = phend[tmp[0]]
        if not tmp[0] in catd:
            catd[tmp[0]] = []
        maj_cats = tmp[1].split(';')
        tmp1 = []
        for i in maj_cats:
            tmp2 = []
            for j in i.split('.'):
                if not j in ["General","Other","Abnormality","Abnormalities","Present"]:
                    tmp2.append(j)
                cat = '.'.join(tmp2)
                tmp1.append(cat)
        catd[tmp[0]] = list(set(tmp1))
catf.close()

# read in gender phenotypes file
genf = open(sys.argv[4],'r')
gend = dict()
for line in genf:
    tmp = line.strip().split()
    gend[tmp[0]] = tmp[1]
genf.close()

# read in mutations file
casesf = open(sys.argv[1],'r')
idx = dict()
mutd = dict()
anno = dict()
for line in casesf:
    tmp = line.strip().split(",")
    if tmp[0] == "Subject": # header line
        for i in xrange(len(tmp)):
            idx[tmp[i]] = i
    else: # data lines
        if not tmp[0] in mutd:
            mutd[tmp[0]] = [[],[],[],[]] # [LGD, D.Mis3, ND.Mis3, other]
        if not tmp[0] in anno:
            anno[tmp[0]] = []
        ec = phend[tmp[0]][0]
        ndd = phend[tmp[0]][1]
        mutgene = tmp[5]
        # get extracardiac major categories for patient if there are any
        if tmp[0] in catd:
            mc = catd[tmp[0]]
        else:
            mc = []  
        type = ''
        # LoF cases
        if tmp[idx["ExonicFunc.refgene"]] in ['stopgain','stoploss','frameshift insertion','frameshift deletion'] or 'splicing' in tmp[idx["Func.refgene"]]: 
            type = "LoF"
            mutd[tmp[0]][0].append(mutgene)
        elif tmp[idx["ExonicFunc.refgene"]] in ['nonsynonymous SNV']:
            # D.Mis3 cases
            if tmp[idx["MetaSVM_pred"]] == "D":
                type = "D.Mis3"
                mutd[tmp[0]][1].append(mutgene)
            # ND.Mis3 cases
            else:
                type = "ND.Mis3"
                mutd[tmp[0]][2].append(mutgene)
        # other mutation cases
        else:
            type = "other"
            mutd[tmp[0]][3].append(mutgene)
        anno[tmp[0]] = [mutd[tmp[0]],mc,ec,ndd]
casesf.close()

cat = str(sys.argv[5])
# calculate summary statistics
lof_dmis_tot = []
counts = [[],[],[]]
ndd_m = 0
ndd_f = 0
non_m = 0
non_f = 0
unc_m = 0
unc_f = 0
for pt in anno:
    nddstat = anno[pt][3]
    ecstat = anno[pt][2]
    mclist = anno[pt][1]
    mutgenes = anno[pt][0][0] + anno[pt][0][1] + anno[pt][0][2] + anno[pt][0][3]
    lofdmisgenes = anno[pt][0][0] + anno[pt][0][1]
    if len(lofdmisgenes) > 0: # if patient has LoF or D.Mis3 mutation
        #if not nddstat == -1: # exclude uncertain cases
        #    lof_dmis_tot.append(pt)
        if len(list(set(lofdmisgenes).intersection(genes))) > 0: # if patient has LoF/D.Mis3 in risk genes list
            # if patient has extracardiac major category
            #if cat in mclist:
            if nddstat == 1:
                counts[0].append((pt, gend[pt]))
                if gend[pt] == "M": ndd_m += 1
                elif gend[pt] == "F": ndd_f += 1
            elif nddstat == 0:
                counts[1].append((pt, gend[pt]))
                if gend[pt] == "M": non_m += 1
                elif gend[pt] == "F": non_f += 1
            elif nddstat == -1:
                counts[2].append((pt, gend[pt]))
                if gend[pt] == "M": unc_m += 1
                elif gend[pt] == "F": unc_f += 1
# add gender information to patient ids in counts list

print sum([len(l) for l in counts])
print [l for l in counts]
print [len(l) for l in counts]
print [(ndd_m, ndd_f), (non_m, non_f), (unc_m, unc_f)]

