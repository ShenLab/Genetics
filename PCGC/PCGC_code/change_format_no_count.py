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

f1 = open('sample_1-00025.hg19_multianno.rare.rare.csv','Ur')
r = csv.reader(f1)
head2 = r.next()
print head2[5:-1]

head = head1 +['Carrier(GT:AB:AD:DP:GQ:MQ0:PL)']+ head2[5:-1] +head3
f = open('1sample_1-00025.csv','wb')
w = csv.writer(f)
w.writerow(head)

for line in r:
    info_all = line[-1].split()
    genotype = [info_all[-3][:3],info_all[-2][:3],info_all[-1][:3]]
    #genotype = info_all[-3:]
    x =info_all[:7]+[';'.join(genotype)]+line[5:-1]+find_correspond_value(info_all[7],head3)

    w.writerow(x)

f1.close()
f.close()


