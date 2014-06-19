def check_rare(l):# find variant with propability < 1% in esp
    for e in l:
        if 'esp.MAF' in e:
            
            freq = float(e.split(',')[-1])
            if freq > 1:                                                                                                       
                return False
            else:
                return True
    return True

def check_GQ(l):
    if './.' in l[0] :
        return False
    else:
        proband= l[0].split(':')
        if len(proband) <=4 :
            return False
        
        proband_GT,proband_GQ = proband[0],int(proband[4])
        if proband_GT == '0/0':
            return False
        if proband_GT == '0/1':
            if l[1] == './.':
                if './.' in l[2] :
                    return proband_GQ > 20
                else:
                    mother_GT = l[2].split(':')[0]
                    if mother_GT != '0/0':
                        return True
                    else:
                        return False
            else:
                father_GT = l[1].split(':')[0]
                if father_GT != '0/0':  #father is 0/1 or 1/1
                    return True
                else:               #father is 0/0, mother needs at least 0/1
                    if './.' in l[2]:
                        return False
                    else:
                        mother_GT = l[2].split(':')[0]
                        if mother_GT == '0/0':
                            return False
                        else:
                            return True
        if proband_GT == '1/1':
            if l[1] == './.':
                if './.' in l[2] :
                    return proband_GQ > 50
                else:
                    mother_GT = l[2].split(':')[0]
                    if mother_GT != '0/0':
                        return proband_GQ > 20                    
                    else:
                        return False
            else:
                father_GT = l[1].split(':')[0]
                if father_GT != '0/0':
                    if './.' in l[2]:
                        return proband_GQ > 20
                    else:
                        mother_GT = l[2].split(':')[0]
                        if mother_GT == '0/0':
                            return False
                        else:
                            return True
                else:
                    return False


def check_deleterious(l):
    for e in l:
        if 'SNPEFF_IMPACT' in e:
            
            impact = e.split('=')[-1]
            if impact == 'LOW':                                                                                                       
                return False
            
    return True




#read  data
f = open('/ifs/scratch/c2b2/ys_lab/yshen/PCGC/analysis/Harvard_Calls/LoF/exome_batch2.vcf')

head = []
for line in f:
    if line[0]=='#': # get  and write head
        if line[1]=='#':
            head.append(line)              
        else:
            head_f = line[1:].split()[:9]
            index = line[1:].split()[9:]
	    #temp = list(index)
            n = len(index)/3
	    add = []
	    for i in range(n):
		if len(index[3*i]) > 7:
			index = index[:i*3]+[index[3*i+1][:-3]]+index[i*3+1:]
			add.extend([3*i])
		if index[3*i]+'-01'!=index[3*i+1]:
			if index[3*i]+'-02' != index[3*i+1]:
				index = index[:i*3+1]+[index[3*i]+'-01']+[index[3*i]+'-02']+index[i*3+1:]
				add.extend([3*i+1,3*i+2])
			else:
				index = index[:i*3+1]+[index[3*i]+'-01']+index[i*3+1:]
				add.extend([3*i+1]) 
		else:
			if index[3*i]+'-02' != index[3*i+2]:
				index = index[:i*3+2]+[index[3*i]+'-02']+index[i*3+2:]
				add.extend([3*i+2])
	    n = len(index)/3
	    '''
            for e in add:
		if e%3 !=0 :
			temp = temp[:e]+['aaa']+temp[e:]
	    for i in range(n):
		print temp[3*i:3*i+3]
	    '''
	    f = [open("temp2/sample_%s.vcf" % index[i*3], "w") for i in range(n)]

            for i in range(n):
                f[i].write("".join(head))
                f[i].write('#'+"\t".join(head_f))
                f[i].write('\t'+index[3*i]+'\t'+index[3*i+1]+'\t'+index[3*i+2]+'\n')
            
    else: #write rare-deleterious-inherited mutations 
        data =  line.split()
        info = data[7].replace(';',' ').split()
	for e in add:
		if e%3 !=0:
			data = data[:e+9]+['./.']+data[e+9:]
        if check_rare(info) and check_deleterious(info):
            for i in range(n):
                temp = data[9+i*3:12+i*3]

                if check_GQ(temp):
                    f[i].write("\t".join(data[:9]+temp)+'\n')                   

for fh in f:
    fh.close()
