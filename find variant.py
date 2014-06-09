import csv
from check_online import check_online_database

print check_online_database('11','48346931')

def find_rare(d):
    # find variant with propability < 1% in 1000 Genomes and GO-ESP
    for e in d:
        if '1KG.score' in e:
            KG_score = float(e[10:])
            if KG_score > 0.01:
                return False
        if 'ESP5400.score' in e :
            ESP5400_score = float(e[14:])
            if ESP5400_score > 0.01:
                return False
    return True 

def check_quality(str1,str2):
    quality1 = str1.split(':')
    quality2 = str2.split(':')
    if len(quality1)<3 or len(quality2)<3 : #no genotype 
        return True
    if float(quality1[3]) > 50 and float(quality2[3]) > 50:# filter those below 50
        return True
    else:
        return False

def check_heterozygous(data_list):
    father_heterozygous = 0
    mother_heterozygous = 0
    if data_list == [[]]:
        return False
    for d in data_list:
        gene1 = d[-2].split(':')[0]
        gene2 = d[-1].split(':')[0]
        if gene1 != '0/0':
            father_heterozygous = 1
        if gene2 != '0/0':
            mother_heterozygous = 1
    if father_heterozygous and mother_heterozygous:
        return True
    else:
        return False
        
    


def write_to_spreadsheet(r,head,info_head):
    temp = []
    
    name = r[-3].split(':')
    sample1 = r[-2].split(':')
    sample2 = r[-1].split(':')
    if r[-2] == './.':
        sample1 = ['./.',' ','  ',' ',' ']
    if r[-1] == './.':
        sample2 = ['./.',' ','  ',' ',' ']
    for i in range(len(name)):
        temp.append(name[i]+'='+sample1[i]+';'+sample2[i])
    r=r[:-3]+ temp  # rearange GT:AD:DP:GQ:PL

    #output to csv for one entry
    w = r[:len(head)]
    for e in info_head:
        j=0
        for i in r:
            if e+'=' in i:
                j=1
                if e == 'DP' and len(i) <= 6:
                    pass
                else:
                    w.append(i[len(e)+1:])
        if j ==0:
            w.append(' ')
                
          
    return w



#read  data
file = open('raw.variants.vcf')

result = []
info_head = []
non_sense = ['functionalClass=stopgainSNV', 'functionalClass=stoplossSNV', 'functionalClass=frameshiftdeletion', 'functionalClass=frameshiftinsertion']
add = 0
SNV=[[]]
gene_previous = ''

for line in file:
    if line[0]=='#': # get head
        if line[1]=='#':
            if line[:6] == '##INFO' or line[:8] == '##FORMAT':
                info_head.append(line.replace(',','=').split('=')[2])               
        else:
            head = line[1:].split()[:7]
    else:
        data =  line.replace(';',' ').split()
        if data[0] in ['X','Y']:# filter sex chromosome 
            break
        if find_rare(data): #find rare
            if check_quality(data[-2],data[-1]):#ensure quality > 50
                for e in non_sense: #find non-sense 
                    if e in line:
                        online = check_online_database(data[0],data[1])
                        print online
                        if not online:
                            add = 1
                if 'functionalClass=nonsynonymousSNV' in line: #find missense, union PolyPhen-2 and SIFT
                    if 'ljb_sift.pred=D' in line or 'ljb_pp2.pred=P' in line or 'ljb_pp2.pred=D'in line:
                        add = 1
                
                if add == 1:
                    for e in data:
                        if 'geneName' in e: #get gene name
                            gene_Name = e[9:]                    
                    if gene_Name == gene_previous:#find those with same gene
                        SNV.append(data)
                        add = 0
                    else:
                        add = 0
                        gene_previous = gene_Name
                        if check_heterozygous(SNV):
                            result.extend(SNV)
                        SNV=[data]
file.close()
if check_heterozygous(SNV):
    result.extend(SNV)

print len(result)     

with open('rare_variant.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(head+info_head)
    for e in result:
          writer.writerow(write_to_spreadsheet(e,head,info_head))

csvfile.close()

