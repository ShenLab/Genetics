## Usage: python parsetableanno2.py eval_domains_training_tabled.hg19_multianno.csv > eval_domains_training_multianno.parsed.csv
## Purpose: parses output of table_annovar.pl, pulls scores from each tool, pulls true labels
## Output: For each variant, for each tool, outputs score
##		format: chr,pos,ref,alt,comments,sift,pp2hdiv,pp2hvar,radialsvm,lr,cadd_phred
import sys

annof = open(sys.argv[1],'r')

colidx = dict()
# header line
print "#chr,pos,ref,alt,comment,sift,pp2hdiv,pp2hvar,radialsvm,lr,cadd_phred,label"
for line in annof:
	pred = [0,0,0,0,0,0] # sift,pp2hdiv,pp2hvar,radialsvm,lr,cadd_phred
	truth = 0
	tmp = line.strip().split(',')
	if tmp[0] == "Chr":
		for i in xrange(len(tmp)):
			colidx[tmp[i]] = i
	else:
		# pull scores
		if not tmp[colidx['SIFT_score']] == '.':
			pred[0] = 1 - float(tmp[colidx['SIFT_score']])
		else:
			pred[0] = tmp[colidx['SIFT_score']]
		pred[1] = tmp[colidx['Polyphen2_HDIV_score']]
		pred[2] = tmp[colidx['Polyphen2_HVAR_score']]
		pred[3] = tmp[colidx['RadialSVM_score']]
		pred[4] = tmp[colidx['LR_score']]
		pred[5] = tmp[colidx['CADD_phred']]
		# fetch true label
		if tmp[colidx['Otherinfo']][1:3] == "TP": truth = 1
		# fetch comment info
		if "|" in tmp[colidx['Otherinfo']]:
			comm = tmp[colidx['Otherinfo']].split("|")[1]
		else:
			comm = tmp[colidx['Otherinfo']].split("\t")[1][:-1]
		# write out chr,pos,ref,alt,comment,sift,pp2hdiv,pp2hvar,radialsvm,lr
		out = [tmp[0],tmp[1],tmp[3],tmp[4],comm] + [str(i) for i in pred] + [str(truth)]
		print ','.join(out)
annof.close()
