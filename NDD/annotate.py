## Usage: python annotate.py cases_metaSVM.csv phenotypecombined.csv kegg_notch.txt kegg_tgfbeta.txt kegg_wnt.txt hmgs.txt txfactor.txt hhe.25.txt hbe.25.txt asd.txt
## Purpose: Creates table of annotations with the columns below. Each row is a subject.  Each column is the count of that mutation type.
## Output: subj ID, LoF mutation info, D.Mis3 mutation info, ND.Mis3 mutation info, indel mutation info, extracardiac, diagnosis
import sys

glists = dict()
## Function to read in list of genes from file and return list object
def parsegenes(fname):
    glist = []
    tmpf = open(fname,'r')
    for line in tmpf:
        tmp = line.strip().split()
        glist.append(tmp[0])
    tmpf.close()
    return glist

## Function to scan gene lists for presence/absence of a given gene and return list object
def scangenes(gldict, gene):
    out = [0]*len(gldict)
    sort = sorted(gldict, key=gldict.get)
    for i in xrange(len(sort)):
        if gene in glists[sort[i]]:
            out[i] = 1
    return out # order of annotations: 'notch', 'wnt', 'hbe', 'txf', 'tgfb', 'hhe', 'asd', 'hmg'
    # order sorted: 'hbe', 'hhe', 'txf', 'tgfb', 'notch', 'hmg', 'wnt', 'asd'

# read gene lists and store in dict as {listname : [genes]}
glists["notch"] = parsegenes(sys.argv[3])
glists["tgfb"] = parsegenes(sys.argv[4])
glists["wnt"] = parsegenes(sys.argv[5])
glists["hmg"] = parsegenes(sys.argv[6])
glists["txf"] = parsegenes(sys.argv[7])
glists["hhe"] = parsegenes(sys.argv[8])
glists["hbe"] = parsegenes(sys.argv[9])
glists["asd"] = parsegenes(sys.argv[10])

# read phenotypes file and store in dict as {id : [ndd, ec]}
phend = dict()
phenf = open(sys.argv[2],'r')
for line in phenf:
    if not line.startswith("ID"):
        tmp = line.strip().split(",")
        id = tmp[0]
        # read in NDD status [NDD, non-NDD, uncertain] => [1, 0, -1]
        if tmp[3] == "NDD": ndd = 1
        elif tmp[3] == "non-NDD": ndd = 0
        elif tmp[3] == "uncertain": ndd = -1
        # read in extracardiac status [Extra-pos, Extra-neg] => [1, 0]
        if tmp[4] == "Extra-pos": ec = 1
        elif tmp[4] == "Extra-neg": ec = 0
        phend[id] = [ndd, ec]
phenf.close()

# read cases (mutations) file and generate annotations
annos = dict() # dict in form: {id : [[mutations]]}
idx = dict() # dictionary to map column names to indexes
f = open(sys.argv[1],'r')
for line in f:
    tmp = line.strip().split(",")
    type = ''
    gene = ''
    if tmp[0] == "Subject":
        for i in xrange(len(tmp)):
            idx[tmp[i]] = i # idx maps column names to indexes
    else:
        if not tmp[0] in annos:
            annos[tmp[0]] = []
        # record mutation information as: [type, gene, 'notch', 'wnt', 'hbe', 'txf', 'tgfb', 'hhe', 'asd', 'hmg']
        # LoF mutations
        if tmp[idx["ExonicFunc.refgene"]] in ['stopgain', 'stoploss', 'frameshift insertion', 'frameshift deletion'] or 'splicing' in tmp[idx["Func.refgene"]]:
            type = "LoF"
        # D.Mis3 and ND.Mis3 mutations
        elif "nonsynonymous SNV" in tmp[idx["ExonicFunc.refgene"]]:
            if tmp[idx["MetaSVM_pred"]] == "D":
                type = "D.Mis3"
            else:
                type = "ND.Mis3"
        # InDel mutations
        elif tmp[idx["ExonicFunc.refgene"]] in ["nonframeshift insertion", "nonframeshift deletion"]:
            type = "indel"
        gene = tmp[idx["Gene"]]
        scanout = scangenes(glists, gene)
        tmpout = [type, gene] + scanout
        annos[tmp[0]].append(tmpout)
f.close()

#print "Subject,LoF_notch,LoF_tgfb,LoF_wnt,LoF_hmg,LoF_txf,LoF_hhe,LoF_hbe,LoF_asd,DMS_notch,DMS_tgfb,DMS_wnt,DMS_hmg,DMS_txf,DMS_hhe,DMS_hbe,DMS_asd,NDMS_notch,NDMS_tgfb,NDMS_wnt,NDMS_hmg,NDMS_txf,NDMS_hhe,NDMS_hbe,NDMS_asd,ID_notch,ID_tgfb,ID_wnt,ID_hmg,ID_txf,ID_hhe,ID_hbe,ID_asd,extracardiac,diagnosis"
print "Subject,LoF_notch,LoF_tgfb,LoF_wnt,LoF_hmg,LoF_txf,LoF_hhe,LoF_hbe,LoF_asd,DMS_notch,DMS_tgfb,DMS_wnt,DMS_hmg,DMS_txf,DMS_hhe,DMS_hbe,DMS_asd,extracardiac,diagnosis"
for subj in annos:
    # features
    LoF_notch, LoF_tgfb, LoF_wnt, LoF_hmg, LoF_txf, LoF_hhe, LoF_hbe, LoF_asd = 0,0,0,0,0,0,0,0
    DMS_notch, DMS_tgfb, DMS_wnt, DMS_hmg, DMS_txf, DMS_hhe, DMS_hbe, DMS_asd = 0,0,0,0,0,0,0,0
    #NDMS_notch, NDMS_tgfb, NDMS_wnt, NDMS_hmg, NDMS_txf, NDMS_hhe, NDMS_hbe, NDMS_asd = 0,0,0,0,0,0,0,0
    #ID_notch, ID_tgfb, ID_wnt, ID_hmg, ID_txf, ID_hhe, ID_hbe, ID_asd = 0,0,0,0,0,0,0,0
    ec = phend[subj][1]
    # label
    diag = phend[subj][0]
    if len(annos[subj]) > 0: # for subjects with mutations
        for mut in annos[subj]:
            # order sorted: 'hbe', 'hhe', 'txf', 'tgfb', 'notch', 'hmg', 'wnt', 'asd'
            if mut[0] == "LoF":
                if mut[2] == 1: LoF_hbe += 1
                if mut[3] == 1: LoF_hhe += 1
                if mut[4] == 1: LoF_txf += 1
                if mut[5] == 1: LoF_tgfb += 1
                if mut[6] == 1: LoF_notch += 1
                if mut[7] == 1: LoF_hmg += 1
                if mut[8] == 1: LoF_wnt += 1
                if mut[9] == 1: LoF_asd += 1
            elif mut[0] == "D.Mis3":
                if mut[2] == 1: DMS_hbe += 1
                if mut[3] == 1: DMS_hhe += 1
                if mut[4] == 1: DMS_txf += 1
                if mut[5] == 1: DMS_tgfb += 1
                if mut[6] == 1: DMS_notch += 1
                if mut[7] == 1: DMS_hmg += 1
                if mut[8] == 1: DMS_wnt += 1
                if mut[9] == 1: DMS_asd += 1
            '''
            elif mut[0] == "ND.Mis3":
                if mut[2] == 1: NDMS_hbe += 1
                if mut[3] == 1: NDMS_hhe += 1
                if mut[4] == 1: NDMS_txf += 1
                if mut[5] == 1: NDMS_tgfb += 1
                if mut[6] == 1: NDMS_notch += 1
                if mut[7] == 1: NDMS_hmg += 1
                if mut[8] == 1: NDMS_wnt += 1
                if mut[9] == 1: NDMS_asd += 1
            elif mut[0] == "indel":
                if mut[2] == 1: ID_hbe += 1
                if mut[3] == 1: ID_hhe += 1
                if mut[4] == 1: ID_txf += 1
                if mut[5] == 1: ID_tgfb += 1
                if mut[6] == 1: ID_notch += 1
                if mut[7] == 1: ID_hmg += 1
                if mut[8] == 1: ID_wnt += 1
                if mut[9] == 1: ID_asd += 1
            '''
    outanno = subj + ',' + ','.join(map(str,[LoF_notch, LoF_tgfb, LoF_wnt, LoF_hmg, LoF_txf, LoF_hhe, LoF_hbe, LoF_asd])) + ',' + ','.join(map(str,[DMS_notch, DMS_tgfb, DMS_wnt, DMS_hmg, DMS_txf, DMS_hhe, DMS_hbe, DMS_asd])) #+ ',' + ','.join(map(str,[NDMS_notch, NDMS_tgfb, NDMS_wnt, NDMS_hmg, NDMS_txf, NDMS_hhe, NDMS_hbe, NDMS_asd])) + ',' + ','.join(map(str,[ID_notch, ID_tgfb, ID_wnt, ID_hmg, ID_txf, ID_hhe, ID_hbe, ID_asd]))
    print outanno + ',' + str(ec) + ',' + str(diag)
