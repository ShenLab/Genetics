def check_rare(l):# find variant with propability < 0.5% in esp
    for e in l:
        if '1KGfreq=' in e:
            freq = float(e.split('=')[-1].split(',')[0])
            if freq > 0.001:
                return False
        if 'ESPfreq=' in e:
	    freq = float(e.split('=')[-1].split(',')[0])
	    if freq > 0.001:
		return False

        if 'AC' in e.split('='):
	    freq = float(e.split('=')[-1].split(',')[0])
	    if freq > 3:
		return False
	if 'AF' in e.split('='):
	    freq = float(e.split('=')[-1].split(',')[0])
	    if freq > 0.015:
		return False
    return True

def check_denovo(l,indel):
    proband= l[0].split(':')
    Father = l[1].split(':')
    Mother = l[2].split(':')
    proband_GT = proband[0]
    Father_GT = Father[0]
    Mother_GT = Mother[0]
    
    if proband_GT in ['./.','0/0']: #proband genotype missing or not denovo
        return False
    else:#GT:AD:DP:GQ:PL	0/1:16,9:.:99:205,0,501
        
        if Father_GT == '0/0' and Mother_GT == '0/0':
            print proband
            proband_GQ = max(map(int,proband[3].split(',')))
            proband_AD = int(proband[1].split(',')[1])
            Father_GQ = max(map(int,Father[3].split(',')))
            if Father[2] =='.':
                Father_DP = sum(map(int,Father[1].split(',')))

            else:
                Father_DP = int(Father[2])

            Mother_GQ = max(map(int,Mother[3].split(',')))
            if Mother[2] =='.':
                Mother_DP = sum(map(int,Mother[1].split(',')))

            else:
                Mother_DP = int(Mother[2])

            if indel:
                return proband_GQ >= 80 and Father_GQ >= 30 and Mother_GQ >= 30 and proband_AD >= 8 and Father_DP >=12 and Mother_DP >= 12
            else:
                return proband_GQ >= 70 and Father_GQ >= 30 and Mother_GQ >= 30 and proband_AD >= 6 and Father_DP >=12 and Mother_DP >= 12
        else:
            return False

f = open('sscpedree')
FORMAT={}
for line in f:
    alll = line.split()
    for e in alll:       
        if e[-2:] == 's1':
            if e[:-2]+'fa' in alll and e[:-2]+'mo' in alll:
                FORMAT[e]=[e[:-2]+'fa',e[:-2]+'mo',alll.index(e),alll.index(e[:-2]+'fa'),alll.index(e[:-2]+'mo')]
f.close()

f = open('SSC_state1_state2.annotated.hardfiltered.vcf')


head = []
for line in f:
    if line[0]=='#': # get and write head
        if line[1]=='#':
            head.append(line)
        else: #
            sample = line[1:].split()
            head_f = sample[:9] # CHROM ... FORMAT
            
            
            f = [open("SSC_state1_state2/"+key+'_'+FORMAT[key][0]+'_'+FORMAT[key][1]+'.vcf', "w") for key in FORMAT]
            n = len(f)
            start=len("SSC_state1_state2/")
            for i in range(n):
                f[i].write("".join(head))
                f[i].write('#'+"\t".join(head_f)+'\t')
                key = f[i].name[start:start+7]
                sscc_index, sscf_index, sscm_index = FORMAT[key][-3:]
                temp = [sample[sscc_index],sample[sscf_index],sample[sscm_index]]
                f[i].write('\t'.join(temp)+'\n')
    
    
    else: 
        data =  line.split()
        indel = not (len(data[4].split(',')[0]) == len(data[3]))
        info = data[7].replace(';',' ').split() # AC=... AF= ...
        if True:
            for i in range(n):
                key = f[i].name[start:start+7]
                sscc_index, sscf_index, sscm_index = FORMAT[key][-3:]
                temp = [data[sscc_index],data[sscf_index],data[sscm_index]]
                if True:
                    f[i].write("\t".join(data[:9]+temp)+'\n')

for fh in f:
    fh.close()

