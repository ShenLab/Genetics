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
            key =  '!'.join(line[:22])
            if key not in mutation:
                mutation[key]=[name]
            else:
                mutation[key].append(name)
        f1.close()
        n += 1
print n

f = open('PCGC_summary.csv','wb')
w = csv.writer(f)
w.writerow(['count','freq']+head)
for key in mutation:
    content= key.split('!')
    count = len(mutation[key])
    sample = ';'.join(mutation[key])
    content.insert(0,float(count)/n)
    content.insert(0,count)
    content.append(sample)
    w.writerow(content)
f.close()

