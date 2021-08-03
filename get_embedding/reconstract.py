#! -*-encoding:utf-8-*-
import numpy as np 
from collections import defaultdict
def gettriples(target,slicen,id2gene,id2disease,slicetype="LOF",switch=None):
    allitem = {**id2disease,**id2gene}
    maxdisease = max(id2disease)
    for row,col in zip(slicen[0],slicen[1]):
        if switch==None:
            function = target[row][col]
            if float(function)>0.53:
                line = "{}\t{}\t{}\t{}\n".format(allitem[col],allitem[row],slicetype,function)
                yield line   
        else:
            function=switch
            line = "{}\t{}\t{}\n".format(allitem[col],allitem[row],function)
            yield line
        
def dumpX(Xpath,outputpath,id2gene,id2disease):
    wf = open(outputpath,'w')
    x = np.load(Xpath)
    lof = x[0].nonzero()
    gof = x[1].nonzero()
    poollof = set()
    for line in gettriples(x,lof,id2gene,id2disease,switch="LOF"):
        wf.write(line)
    for line in gettriples(x,gof,id2gene,id2disease,switch="GOF"):
        wf.write(line)

def dumpARAT(Apath,Rpath,outputpath,id2gene,id2disease):
    wf = open(outputpath,'w')
    A = np.load(Apath)
    R = np.load(Rpath)
    ARAT = [np.dot(np.dot(A,R[0]),A.T),np.dot(np.dot(A,R[1]),A.T)]
    ARAT_LOF = ARAT[0].nonzero()
    ARAT_GOF = ARAT[1].nonzero()
    for line in gettriples(ARAT[0],ARAT_LOF,id2gene,id2disease,slicetype="LOF",switch=None):
        wf.write(line)
    for line in gettriples(ARAT[1],ARAT_GOF,id2gene,id2disease,slicetype="GOF",switch=None):
        wf.write(line)

def dumptrain_or_total(inputfile):
    rf = open(inputfile,'r')
    hashdir = defaultdict(list)
    for line in rf:
        contents = line.strip().split('\t')
        genes = contents[0].split(';')
        for g in genes:
            func = contents[1]
            mesh = "MESH:"+contents[2]
            tuple_ = tuple(sorted([g,mesh]))
            hashdir[tuple_].append(func)
    for t,l in hashdir.items():
        if len(l)>=2:
            hashdir[t]="COM"
        else:
            hashdir[t]=l[0]
    return hashdir            

def uniq(inputfile):
    hashdir = defaultdict(list)
    rf = open(inputfile,"r")
    for line in rf:
        contents = line.strip().split("\t")
        item1 = contents[0]
        item2 = contents[1]
        tuple_ = tuple(sorted([item1,item2]))
        hashdir[tuple_].append(contents[2])
    for t,l in hashdir.items():
        if len(l)==2:
            hashdir[t]="COM"
        else:
            hashdir[t]=l[0]
    return hashdir

def mapping(path):
    hashdir = {}
    with open(path) as rf:
        for line in rf:
            contents = line.strip().split("\t")
            id_ = contents[1]
            name = contents[0]
            hashdir[id_]=name
    return hashdir

def main():
    Xpath = "output/X.npy"
    Apath = "output/A.npy"
    Rpath = "output/R.npy"
    Xoutputpath = "analysis/X.txt"
    ARAToutputpath = "analysis/ARAT.txt"
    gene2id = np.load("output/gene2index.npy").item()
    disease2id = np.load("output/disease2index.npy").item()
    id2gene = {id:gene for gene,id in gene2id.items()}
    id2disease = {id:disease for disease,id in disease2id.items()}
    dumpX(Xpath,Xoutputpath,id2gene,id2disease)
    dumpARAT(Apath,Rpath,ARAToutputpath,id2gene,id2disease)
    newxdir = uniq("analysis/X.txt")
    newARATdir = uniq("analysis/ARAT.txt")
    traindir = dumptrain_or_total("../DataForWA/gene_function_disease_total_train.txt")
    totaldir = dumptrain_or_total("../DataForWA/gene_function_disease_total_done.txt")
    count=0
    predicts={}
    for gene_mesh,func in newARATdir.items():

        xfunc = newxdir.get(gene_mesh,-1)
        if gene_mesh in newxdir and xfunc==func:
            count+=1
        else:
            predicts[gene_mesh]=func

    mesh2symbol=mapping("../DataForWA/diseasename2mesh.txt",)
    genid2symbol = mapping("../DataForWA/genename2geneid.txt")
    pref = open("predicts.txt",'w')
    for key,value in predicts.items():
        geneid = key[0]
        if geneid.startswith("MESH"):
            pass
        else:
            disease = key[1].split(":")[1]
            function = value
            # print(geneid,disease,function)
            pref.write("{}\t{}\t{}\n".format(genid2symbol[geneid],function,mesh2symbol[disease]))
        

    print("Train Recall:{}".format(count/len(newxdir)))
    print("Train Presion:{}".format(count/len(newARATdir)))

    testdir = {key:value for key,value in totaldir.items()-traindir.items()}

    count2 = 0
    # print(predicts)
    for gene_mesh,func in predicts.items():
        xfunc = testdir.get(gene_mesh,-1)
        if gene_mesh in testdir and xfunc==func:
            count2+=1

    print("Test Recall:{}".format(count2/len(testdir)))
    print("Test Presion:{}".format(count2/len(predicts)))
    print("X:{}".format(len(newxdir)))
    print("ARAT:{}".format(len(newARATdir)))
    print("Test:{}".format(len(testdir)))
    print("Predict:{}".format(len(predicts)))
if __name__=="__main__":
    main()
