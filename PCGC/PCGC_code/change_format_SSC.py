import csv


def find_correspond_value(s,head):
    info = s.split(';')
    value_l = []

    
    for h in head:
        add = 0
        for i in info:
            entry = i.split('=')
            if len(entry) !=1 and entry[0]==h:
                add = 1
                
                value_l.append(entry[-1])

        if add == 0:
            value_l.append('.')
    return value_l
        

#read  data

f3 = open('head.txt','r')

head1 = f3.readline().split()
head3 = f3.readline().split()

f3.close()

f1 = open('recessive.csv','Ur')
r = csv.reader(f1)
head2 = r.next()
head = ['count', 'freq'] + head1 +['Carrier(GT)']+['Inherited from ']+['Carrier AD:GQ']+ head2[7:-1] +head3
f = open('SSC_recessive_readable.csv','wb')
w = csv.writer(f)
w.writerow(head)


proband_list =[]
parents_list=[]
proband_AD_list=[]

for line in r:
    info_all = line[-1].split()
    
    for e in line[-2].split(');'):

        ID, info =  e.split('(')
        proband_ID,father_ID,mother_ID = ID.split('_')
        sample_info = info[:-1].split(';')

        proband_geno = sample_info[0][:3]
        Father_geno = sample_info[1][:3]
        Mother_geno = sample_info[2][:3]
        genoinfo = sample_info[0].split(':')

        if len(genoinfo) == 7:  #AD:GQ
            proband_geno_AD = genoinfo[2]+':'+genoinfo[-3]
        else:
            proband_geno_AD = genoinfo[1]+':'+genoinfo[-3]
            
        if proband_geno == '1/1':
            parents_list.append(father_ID+'('+Father_geno+')'+','+mother_ID+'('+Mother_geno+')')
        if proband_geno == '0/1':
            if '1' in Father_geno: # 0/1 0/1 x/x
                parents_list.append(father_ID+'('+Father_geno+')')
            elif '0/0' in Father_geno: # 0/1 0/0 0/1 or 1/1 or ./.
                parents_list.append(mother_ID+'('+Mother_geno+')')
            else:
                if '1' in Mother_geno:# 0/1 ./. 0/1 or 1/1
                    parents_list.append(mother_ID+'('+Mother_geno+')')
                elif '0/0' in Mother_geno:# 0/1 ./. 0/0
                    parents_list.append(father_ID+'('+Father_geno+')')
                else: #0/1 ./. ./.
                    parents_list.append(father_ID +'('+Father_geno+')'+','+mother_ID +'('+Mother_geno+')')

        
        proband_list.append(proband_ID+'('+proband_geno+')')
        proband_AD_list.append(proband_geno_AD)
        #print proband_geno_AD
    geno_info_trio= [';'.join(proband_list)]+[';'.join(parents_list)]+[';'.join(proband_AD_list)]
    x =line[:2]+info_all[:7]+geno_info_trio+line[7:-2]+find_correspond_value(info_all[7],head3)

    w.writerow(x)
    proband_list =[]
    proband_AD_list = []
    parents_list = []

f1.close()
f.close()
