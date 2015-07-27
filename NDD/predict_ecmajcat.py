## Usage: python predict6.py cases_metaSVM.csv phenotypecombined.csv phenotypes_extracardiac_mb_wkc.csv <Gene List> 
## Purpose: Parses patient mutation information, NDD status, extracardiac major categories; 
##            calculates posterior probability of NDD for each major category
## Output: Sensitivity and Counts for NDD, non-NDD, uncertain among patients with damaging mutations in risk genes list 
##            AND who show EC manifestation from list of major categories with posterior > cutoff
import sys

## open files
casesf = open(sys.argv[1],'r')
phenf = open(sys.argv[2],'r')
extraf = open(sys.argv[3],'r')
genef = open(sys.argv[4],'r')

## process risk genes file
tmp = []
for line in genef:
    if not line.startswith("#"):
        tmp.append(line.strip().split(',')[0])
genef.close()
genes = list(set(tmp)) # delete duplicates

## process phenotypes file
phend = dict()
for line in phenf:
    tmp = line.strip().split(",")
    if not tmp[0] == "ID":
        # capture NDD diagnosis
        if tmp[3] == "NDD": ndd = 1
        elif tmp[3] == "non-NDD": ndd = 0
        elif tmp[3] == "uncertain": ndd = -1
        # capture extracardiac status
        if tmp[4] == "Extra-pos": ec = 1
        elif tmp[4] == "Extra-neg": ec = 0
        # record extracardiac and NDD status for each patient
        phend[tmp[0]] = [ec, ndd]
phenf.close()

## process phenotypes extracardiac file
def processextra(extrafile):
	tmp_catd = dict() # dictionary to hold patient IDs for each major category 
	tmp_patd = dict() # dictionary to hold major categories for each patient ID
	for line in extrafile:
		tmp = line.strip().split(',')
		if not tmp[0] == "ID":
			nddstat = phend[tmp[0]][1]
			if not tmp[0] in tmp_patd:
				tmp_patd[tmp[0]] = [[],nddstat]
			maj_cats = tmp[1].split(';')
			tmp1 = []
			for i in maj_cats:
				tmp2 = []
				for j in i.split('.'):
					if not j in ["General","Other", "Abnormality", "Abnormalities", "Present"]:
						tmp2.append(j)
				cat = '.'.join(tmp2)
				if not cat in tmp_catd:
					tmp_catd[cat] = [[],[],[]]
				if nddstat == 1: tmp_catd[cat][0].append(tmp[0])
				elif nddstat == 0: tmp_catd[cat][1].append(tmp[0])
				elif nddstat == -1: tmp_catd[cat][2].append(tmp[0])
				tmp1.append(cat)
			tmp_patd[tmp[0]][0] = list(set(tmp1))
	return(tmp_catd, tmp_patd)
topcats = [] # list of categories with posteriors above the cutoff
outd = dict()
catd, patd = processextra(extraf)
for cat in catd:
		a, b = 1, 1
		stats = [len(list(set(catd[cat][0]))), len(list(set(catd[cat][1]))), len(list(set(catd[cat][2])))]
		posterior = float(stats[0]+a)/float(stats[0]+stats[1]+a+b)
		outd[cat] = posterior
		if posterior >= 0.6: # set cutoff here
			topcats.append(cat) 
extraf.close()

## process mutations file
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
        if tmp[0] in patd:
            mc = patd[tmp[0]][0]
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

## calculate summary statistics
counts = [[],[],[]]
for pt in anno:
    nddstat = anno[pt][3]
    ecstat = anno[pt][2]
    mclist = anno[pt][1]
    lofdmisgenes = anno[pt][0][0] + anno[pt][0][1]
    if len(lofdmisgenes) > 0: # if patient has LoF or D.Mis3 mutation
        if len(list(set(lofdmisgenes).intersection(genes))) > 0: # if patient has LoF/D.Mis3 in risk genes list
            if len(list(set(mclist).intersection(topcats))) > 0: # if patient has extracardiac major cat with posterior >= cutoff
                if nddstat == 1:
                    counts[0].append(pt)
                elif nddstat == 0:
                    counts[1].append(pt)
                elif nddstat == -1:
                    counts[2].append(pt)
print sum([len(l) for l in counts])
print [len(l) for l in counts]
sens = float(len(counts[0]))/float(len(counts[0])+len(counts[1]))
print sens


