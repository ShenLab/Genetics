import csv
import os

def check_rare(line):
    # find variant with propability < 1% in 1000 Genomes and GO-ESP
    esp6500si_all = line[esp_index]
    kg = line[kg_index]
    
    if esp6500si_all == 'NA' :
        if kg == 'NA':
            return True
        else:
            kg = float(kg)
            if kg > 0.01:
                return False
            else:
                return True
    else:
        esp6500si_all = float(esp6500si_all)
        if esp6500si_all > 0.01:
            return False
        else:
            if kg == 'NA':         
                return True
            else:
                kg = float(kg)
                if kg > 0.01:
                    return False
                else:
                    return True        

def check_deleterious(line):
    LJB_sift_pre = line[LJB_sift]
    LJB_pp2_pre = line[LJB_pp2]
    if LJB_sift_pre in ['NA','.'] or LJB_pp2_pre == 'NA':
        
        nonsynonymous_d = False
    else:
        if LJB_sift_pre == 'D' or LJB_pp2_pre in['D','P']:
            nonsynonymous_d = True
        else:
            nonsynonymous_d = False
    if line[ExonicFunc_index] in non_sense or 'splicing' in line[ExonicFunc_index-2] or nonsynonymous_d:
            return True
    else:
        return False


#read  data

for e in  os.listdir('afterannovar'):
    if e.endswith(".csv"):
        f1 = open('afterannovar/'+e,'Ur')
        r = csv.reader(f1)
        f2 = open('PCGC_rare_csv/'+e[:-3]+'rare.csv','wb')
        w = csv.writer(f2)
        head = r.next()
        w.writerow(head)
        index = 0
        for e in head:
            if e == 'esp6500si_all':
                esp_index = index
            if e == '1000g2012apr_all':
                kg_index = index
            if e == 'ExonicFunc.refGene':
                ExonicFunc_index = index
            if e == 'LJB_SIFT_Pred':
                LJB_sift = index
            if e == 'LJB_PolyPhen2_Pred':
                LJB_pp2 = index
            index += 1

        non_sense = ['stopgain SNV', 'stoploss SNV', 'frameshift deletion', 'frameshift insertion']


        for line in r:
            if check_rare(line) and check_deleterious(line) :
                w.writerow(line)
        
            

        f1.close()
        f2.close()
