import csv
import os




#read  data
mutation = {}
n=0
for e in os.listdir('PCGC_rare_csv'):
    if e.endswith(".rare.csv"):
        f1 = open('PCGC_rare_csv/'+e,'Ur')
        r = csv.reader(f1)
        head = r.next()
        name = e[7:14]
        
        for line in r:
            data = line[-1].split()            
            info,genotype =  data[:-4],data[-3:]
            key =  '!'.join(line[:-1])
            sample_geno = name +'('+';'.join(genotype)+')'
            if key not in mutation:
		#mutation[key] = [name]
                mutation[key]=[info,sample_geno]
            else:
		#mutation[key].append(name)
                mutation[key].append(sample_geno)
        f1.close()
        n += 1
print n

f = open('PCGC_summary.csv','wb')
w = csv.writer(f)
w.writerow(['count','freq']+head)
for key in mutation:
    content= key.split('!')
    count = len(mutation[key])-1
    sample = ';'.join(mutation[key][1:])
    content.append(sample)
    content.append('\t'.join(mutation[key][0]))
    content.insert(0,float(count)/n)
    content.insert(0,count)
    w.writerow(content)
f.close()

