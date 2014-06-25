import csv
import os

def check_heterozygous(l):
    father_heterozygous = 0
    mother_heterozygous = 0
    if l == []:
        return False
    for d in l:
        gene1 = d[0]
        gene2 = d[1]
        if gene1 != '0/0':
            father_heterozygous = 1
        if gene2 != '0/0':
            mother_heterozygous = 1
    if father_heterozygous and mother_heterozygous:
        return True
    else:
        return False
    
def check_ambiguous(F,M): #return True if not ambigous like 0/0 0/1
    if '1' in F and '1' in M:
        return False
    else:
        return True
        


for e in  os.listdir('PCGC_rare_csv'):

    if e.endswith(".csv"):

        #open file, write head 
        f1 = open('PCGC_rare_csv/'+e,'Ur')

        r = csv.reader(f1)

        f2 = open('PCGC_recessive_csv/'+e[:-3]+'rare.csv','wb')

        w = csv.writer(f2)

        head = r.next()

        w.writerow(head)

        
        #find gene name index
        index = 0
        for e in head:
            if e =='Gene.refGene':
                gene_name_index = index
                break
            index += 1


        #write homozygous and compound heterozygous 
        gene_previous = ''
        SNV = []
        gene_trio=[]

        for line in r:
            
            gene_name = line[gene_name_index].split('(')[0]
            info_all = line[-1].split()
            proband = info_all[-3].split(':')[0]
            Father = info_all[-2].split(':')[0]
            Mother = info_all[-1].split(':')[0]
            
            if Father == './.':
                if '1' in Mother:
                    Father = '0/0'
                elif Mother == './.':
                    pass
                else:
                    Father = '0/1'
            if Mother == './.':
                if '1' in Father:
                    Mother = '0/0'
                elif Mother == './.':
                    pass
                else:
                    Mother = '0/1'

            if proband == '1/1': # include homozygous 
                w.writerow(line)
            else: #include compound heterozygous
                if check_ambiguous(Father,Mother):    
                    if gene_name == gene_previous:
                        gene_trio.append([Father,Mother])
                        SNV.append(line)
                    else:
                        gene_previous = gene_name
                        if check_heterozygous(gene_trio):
                            for e in SNV:
                                w.writerow(e)
                        SNV = [line]
                        gene_trio = [[Father,Mother]]

            
        if check_heterozygous(gene_trio):
            for e in SNV:
                w.writerow(e)

        f1.close()
        f2.close()

