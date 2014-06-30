def check_rare(l):# find variant with propability < 0.5% in esp
    for e in l:
        if 'esp.MAF' in e:
            freq = float(e.split(',')[-1])
            if freq > 0.5:
                return False
            else:
                return True
    return True

def check_GQ(l):# l =[f,m,s]
    proband = l[0].split(':')
    Father = l[1].split(':')
    Mother = l[2].split(':')
    proband_GT = proband[0]
    Father_GT = Father[0]
    Mother_GT = Mother[0]
    one_missing_one_place = 50
    one_missing_two_place = 30
    two_missing_two_place = 80
    
    
    if './.' in proband_GT : #proband genotype missing
        return False

    else:
        proband_GQ = int(proband[-3])
        
        if proband_GT == '0/0':
            return False
        
        if proband_GT == '0/1':
            if Father_GT == './.':
                if Mother_GT == './.' :
                    return proband_GQ > one_missing_two_place # 0/1 ./. ./.
                else:
                    
                    if Mother_GT != '0/0':# 0/1 ./. 0/1
                        return True
                    else:
                        return proband_GQ > one_missing_one_place # 0/1 ./. 0/0
            else:
                if Father_GT != '0/0':  #0/1  (0/1 or 1/1)
                    return True
                else:               #father is 0/0, mother needs at least 0/1
                    if './.' in Mother_GT: # 0/1 0/0 ./.
                        return proband_GQ > one_missing_one_place
                    else:
                        if Mother_GT == '0/0': # 0/1 0/0 0/0
                            return False
                        else:
                            return True
                                
        if proband_GT == '1/1':
            if Father_GT == './.':
                if Mother_GT == './.' :
                    return proband_GQ > two_missing_two_place # 1/1 ./. ./.
                else:
                
                    if Mother_GT != '0/0': # 1/1 ./. 0/1 or 1/1
                        return proband_GQ > one_missing_one_place
                    else: #1/1 ./. 0/0
                        return False
            else:
                if Father_GT != '0/0':
                    if Mother_GT == './.': # 1/1 0/1 ./.
                        return proband_GQ > one_missing_one_place
                    else:
                        if Mother_GT == '0/0':
                            return False
                        else:
                            return True
                else: # 1/1 0/0 mother
                    return False


def check_deleterious(l):
    for e in l:
        if 'SNPEFF_IMPACT' in e:
            impact = e.split('=')[-1]
            if impact == 'LOW':                                                                                                       
                return False
    return True




#read  data
f = open('/ifs/scratch/c2b2/ys_lab/yshen/PCGC/analysis/Harvard_Calls/denovo/ssc-yale-state1.snpeff.vcf')

head = []
for line in f:
    if line[0]=='#': # get and write head
        if line[1]=='#':
            head.append(line)              
        else: #
            head_f = line[1:].split()[:9] # CHROM ... FORMAT
            index = line[1:].split()[9:] # 11008.fa	11008.mo	11008.s1 ...

            n = len(index)/3
            f = [open("SSC_vcf/" + index[i*3+2]+'_'+index[i*3]+'_'+index[i*3+1] +'.vcf', "w") for i in range(n)]
            for i in range(n):
                f[i].write("".join(head))
                f[i].write('#'+"\t".join(head_f)+'\t')
                f[i].write('\t'.join([index[3*i+2]]+index[3*i:3*i+2])+'\n')
            
            
    else: #write rare-deleterious-inherited mutations 
        data =  line.split()
        info = data[7].replace(';',' ').split()
        if check_rare(info) and check_deleterious(info) and data[6] == 'PASS':
            for i in range(n):
                temp = [data[11+3*i]] + data[9+i*3:11+i*3]

                if check_GQ(temp):
                    f[i].write("\t".join(data[:9]+temp)+'\n')                   

for fh in f:
    fh.close()
 
