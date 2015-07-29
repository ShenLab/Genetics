## Usage: python highrisk.py cases_metaSVM.csv phenotypecombined.csv phenotypes_extracardiac_mb_wkc.csv riskgenes.txt <cutoff> > highriskannos.csv
## Purpose: Generates annotation table containing the following information for each patient:
## 			1) # of damaging (LGD/Deleterious Missense) mutations in genes in riskgenes.txt (default: 0)
## 			2) extracardiac manifestation major category posterior probability of NDD (default: 0, ranges [0,1])
## 			3) diagnosis status (1: NDD, 0: non-NDD, -1: uncertain)
## Output: table of annotations with rows: subjectID, #damaging mutations, #extracardiac posterior, NDD status

import sys

mutf = open(sys.argv[1],'r')
phenf = open(sys.argv[2],'r')
extraf = open(sys.argv[3],'r')
genef = open(sys.argv[4],'r')
cutoff = sys.argv[5]

## process risk genes file
tmp = []
for line in genef:
	tmp.append(line.strip())
genef.close()
genes = list(set(tmp)) # list of unique risk genes

## process phenotypes (diagnosis) file
phend = dict()
for line in phenf:
	tmp = line.strip().split(',')
	if not tmp[0] == "ID":
		# capture NDD status (diagnosis)
		if tmp[3] == "NDD": 
			ndd = 1
		elif tmp[3] == "non-NDD": 
			ndd = 0
		elif tmp[3] == "uncertain": 
			ndd = -1
		phend[tmp[0]] = ndd
phenf.close()

## process extracardiac manifestations file
def processextra(extrafile):
	tmp_catd = dict()
	tmp_patd = dict()
	for line in extrafile:
		tmp = line.strip().split(',')
		if not tmp[0] == "ID":
			ndd = phend[tmp[0]]
			if not tmp[0] in tmp_patd:
				tmp_patd[tmp[0]] = [[], ndd]
			maj_cats = tmp[1].split(';')
			tmp1 = [] # to hold extracardiac major categories for a given patient
			for i in maj_cats:
				tmp2 = [] # to hold extracardiac major categories for processing
				for j in i.split('.'):
					if not j in ["General", "Other", "Abnormality", "Abnormalities", "Present"]:
						tmp2.append(j)
				cat = '.'.join(tmp2)
				if not cat in tmp_catd:
					tmp_catd[cat] = [[],[],[]]
				if ndd == 1: 
					tmp_catd[cat][0].append(tmp[0])
				elif ndd == 0: 
					tmp_catd[cat][1].append(tmp[0])
				elif ndd == -1:
					tmp_catd[cat][2].append(tmp[0])
				tmp1.append(cat)
			tmp_patd[tmp[0]][0] = list(set(tmp1))
	return(tmp_catd, tmp_patd)
# catd holds #NDD/non-NDD/uncertain patients for each extracardiac major category
# patd holds the extracardiac major categories each patient has
catd, patd = processextra(extraf)
## calculate extracardiac major categories posterior probability of NDD
topcats = [] # hold major categories with posterior probabilities >= cutoff
postd = dict() # hold posterior probability values for each major category
for cat in catd:
	a, b = 1, 1
	stats = [len(list(set(catd[cat][0]))), len(list(set(catd[cat][1]))), len(list(set(catd[cat][2])))]
	posterior = float(stats[0]+a)/float(stats[0]+stats[1]+a+b)
	postd[cat] = posterior
	if posterior >= float(cutoff): # set cutoff here
		topcats.append(cat) 
extraf.close()

## process mutations file
idx = dict()
mutd = dict() # for each patient, hold all genes mutated
anno = dict() 
for line in mutf:
	tmp = line.strip().split(',')
	if tmp[0] == "Subject": # handle header line
		for i in xrange(len(tmp)):
			idx[tmp[i]] = i
	else: # data lines
		if not tmp[0] in mutd:
			mutd[tmp[0]] = [[],[],[],[]] # [LGD, Dmis, NDmis, other]
		if not tmp[0] in anno:
			anno[tmp[0]] = []
		ndd = phend[tmp[0]]
		mutgene = tmp[5]
		# get extracardiac major categories for patient if they have any
		if tmp[0] in patd:
			mc = patd[tmp[0]][0] # list of major categories for a given patient
			posterior = 0
			# take highest posterior probability value among all major categories
			for c in mc:
				if postd[c] > posterior:
					posterior = postd[c]
		else:
			mc = []
			posterior = 0
		# LoF mutations
		if tmp[idx["ExonicFunc.refgene"]] in ['stopgain','stoploss','frameshift insertion','frameshift deletion'] or 'splicing' in tmp[idx["Func.refgene"]]: 
			mutd[tmp[0]][0].append(mutgene)
		# Dmis and NDmis mutations
		elif tmp[idx["ExonicFunc.refgene"]] in ['nonsynonymous SNV']:
			if tmp[idx["MetaSVM_pred"]] == "D":
				mutd[tmp[0]][1].append(mutgene)
			else:
				mutd[tmp[0]][2].append(mutgene)
		else:
			mutd[tmp[0]][3].append(mutgene)
		anno[tmp[0]] = [mutd[tmp[0]], mc, posterior, ndd]
mutf.close()

## generate annotations
counts = [0,0,0]
for pt in anno:
	damagingct = 0
	mcposterior = anno[pt][2]
	ndd = anno[pt][3]
	mclist = anno[pt][1]
	lofdmisgenes = anno[pt][0][0] + anno[pt][0][1]
	for g in lofdmisgenes:
		if g in genes:
			damagingct += 1
	print ','.join([pt, str(damagingct), str(mcposterior), str(ndd)])
	'''
	if damagingct > 0 and mcposterior >= float(cutoff):
		print ','.join([pt, str(damagingct), str(mcposterior), str(ndd)])
		if ndd == 1:
			counts[0] += 1
		elif ndd == 0:
			counts[1] += 1
		elif ndd == -1:
			counts[2] += 1
print counts
'''
