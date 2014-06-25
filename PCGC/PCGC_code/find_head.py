import csv
import os




#read  data

f1 = open('/ifs/scratch/c2b2/ys_lab/yshen/PCGC/analysis/Harvard_Calls/LoF/')


head=[]
for line in f1:
    if line[0]=='#':
        if line[:10]=='##INFO=<ID':
            head.append(line.replace(',','=').split('=')[2])
        else:
            head1 = line[1:].split()[:7]
f1.close()

print head
print head1


f = open('head.txt','wb')
f.write('\t'.join(head1)+'\n')
f.write('\t'.join(head))

f.close()

