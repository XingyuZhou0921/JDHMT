#! usr/bin/env python3
# -*- coding:utf-8 -*-
"""
@Author:zhoukaiyin
This version use relation to filter unrelated data.
Here my thoughts are, if gene and disease in triples, it will be our target gene and disease
besides, the gene relate with those gene and those disease will be used.
Don't use CNV,MC
"""
from RescalSideV1 import RescalSide,Rescal,PCA
from scipy.sparse import csr_matrix
import json
import os
import argparse
import numpy as np
from collections import defaultdict 

class Mapping(object):
    def __init__(self,path):
        self.path = path
    
    def convert2id(self):
        gene2id = {}
        r = open(self.path)
        for line in r:
            gene = line.strip().split('\t')[-1]
            gene2id[gene]=len(gene2id)
        return gene2id

    def convert2entity(self):
        gene2id = self.convert2id()
        id2gene = {id:gene for gene,id in gene2id.items()}
        return id2gene

class SideMapping():
    def __init__(self,dictionary,alltriples,filter_flag=False):
        self.dictionary = dictionary
        self.alltriples = alltriples
        self.dc = json.load(open(os.path.join(self.dictionary, "disease_chemical.json")))
        self.dd = json.load(open(os.path.join(self.dictionary, "disease_disease.json")))
        self.dm = json.load(open(os.path.join(self.dictionary, "disease_mutation.json")))
        self.dpa = json.load(open(os.path.join(self.dictionary, "disease_pathway.json")))
        self.dpe = json.load(open(os.path.join(self.dictionary, "disease_phe.json")))
        self.dg = json.load(open(os.path.join(self.dictionary, "disease_gene.json")))

        self.gc = json.load(open(os.path.join(self.dictionary, "gene_chemical.json")))
        self.gg = json.load(open(os.path.join(self.dictionary, "gene_gene.json")))
        self.gm = json.load(open(os.path.join(self.dictionary, "gene_mutation.json")))
        self.gpa = json.load(open(os.path.join(self.dictionary, "gene_pathway.json")))
        self.gpe = json.load(open(os.path.join(self.dictionary, "gene_phe.json")))
        self.gd = json.load(open(os.path.join(self.dictionary, "gene_disease.json")))

        importantdiseases=[]
        importantgenes=[]
        if filter_flag:
            importantdiseases,importantgenes = self.get_important_gene_and_disease(self.alltriples,self.gg,self.dd)
        self.pipline(importantdiseases,importantgenes,filter_flag)

    def get_important_gene_and_disease(self,triples,gg,dd):
        importantdiseases=set()
        importantgenes = set()
        rf = open(self.alltriples)
        for line in rf:
            contents = line.strip().split() 
            disease = contents[2]
            genes = contents[0].split(";")
            for gene in genes:
                importantgenes.add(gene)
                if gene in gg:
                    importantgenes = importantgenes | set(gg[gene])
            importantdiseases.add(disease)
            ndisease = "MESH:"+disease
            if ndisease in dd:
                importantdiseases = importantdiseases | set([d.split(":")[-1] for d in dd[ndisease]])        
        return importantdiseases,importantgenes

    def pipline(self,importantdiseases,importantgenes,filter_flag):
        if filter_flag:
            self.dc = {key:tuple(value) for key,value in self.dc.items() if key.split(":")[-1] in importantdiseases}
            self.dd = {key:tuple(value) for key,value in self.dd.items() if key.split(":")[-1] in importantdiseases}
            self.dm = {key:tuple(value) for key,value in self.dm.items() if key.split(":")[-1] in importantdiseases}
            self.dpa = {key:tuple(value) for key,value in self.dpa.items() if key.split(":")[-1] in importantdiseases}
            self.dpe = {key:tuple(value) for key,value in self.dpe.items() if key.split(":")[-1] in importantdiseases}
            self.dg = {key:tuple(value) for key,value in self.dg.items() if value in importantgenes and key.split(":")[-1] in importantdiseases}

            self.gc = {key:tuple(value) for key,value in self.gc.items() if key in importantgenes}
            self.gg = {key:tuple(value) for key,value in self.gg.items() if key in importantgenes}
            self.gm = {key:tuple(value) for key,value in self.gm.items() if key in importantgenes}
            self.gpa = {key:tuple(value) for key,value in self.gpa.items() if key in importantgenes}
            self.gpe = {key:tuple(value) for key,value in self.gpe.items() if key in importantgenes}
            self.gd = {key:tuple(value) for key,value in self.gd.items() if key in importantgenes and value.split(":")[-1] in importantdiseases}

        else:
            self.dc = {key:tuple(value) for key,value in self.dc.items()}
            self.dd = {key:tuple(value) for key,value in self.dd.items()}
            self.dm = {key:tuple(value) for key,value in self.dm.items()}
            self.dpa = {key:tuple(value) for key,value in self.dpa.items()}
            self.dpe = {key:tuple(value) for key,value in self.dpe.items()}
            self.dg = {key:tuple(value) for key,value in self.dg.items()}


            self.gc = {key:tuple(value) for key,value in self.gc.items()}
            self.gg = {key:tuple(value) for key,value in self.gg.items()}
            self.gm = {key:tuple(value) for key,value in self.gm.items()}
            self.gpa = {key:tuple(value) for key,value in self.gpa.items()}
            self.gpe = {key:tuple(value) for key,value in self.gpe.items()}
            self.gd = {key:tuple(value) for key,value in self.gd.items()}
            
    def hebing(self,diseasdir, genedir):
        entity = []
        for item in diseasdir.values():
            entity.extend(item)
        for item1 in genedir.values():
            entity.extend(item1)
        return sorted(set(entity))

    def merge_all(self,ec,em,epa,epe,ee,e1e2):
        entity = [k for d in [ec,em,epa,epe,ee,e1e2] for k in d.keys()]
        for item2 in ee.values():
            entity.extend(item2)
        return sorted(set(entity))

    def entity2id(self):
        chemicals= self.hebing(self.dc,self.gc)
        mutations = self.hebing(self.dm,self.gm)
        pathways = self.hebing(self.dpa,self.gpa)
        phenotypes = self.hebing(self.dpe,self.gpe)
        diseases = self.merge_all(self.dc,self.dm,self.dpa,self.dpe,self.dd,self.dg)
        genes = self.merge_all(self.gc,self.gm,self.gpa,self.gpe,self.gg,self.gd)
        # map to id
        chemical2id = {item:i for i,item in enumerate(chemicals)}
        mutation2id = {item:i for i,item in enumerate(mutations)}
        pathway2id = {item:i for i,item in enumerate(pathways)}
        phenotype2id = {item:i for i,item in enumerate(phenotypes)}
        gene2id = {item:i for i,item in enumerate(genes)}
        disease2id = {item:i for i,item in enumerate(diseases)}
        return chemical2id,mutation2id,pathway2id,phenotype2id,gene2id,disease2id

    def id2entity(self):
        chemical2id,mutation2id,pathway2id,phenotype2id,gene2id,disease2id = self.entity2id()
        id2chemical={i:c for c,i in enumerate(chemical2id)}
        id2mutation={i:m for m,i in enumerate(mutation2id)}
        id2pathway={i:p for p,i in enumerate(pathway2id)}
        id2phenotype={i:ph for ph,i in enumerate(phenotype2id)}
        id2gene={i:g for g,i in enumerate(gene2id)}
        id2disease={i:d for d,i in enumerate(disease2id)}
        return id2chemical,id2mutation,id2pathway,id2phenotype,id2gene,id2disease

    def counter(self):
        print("Disease-Chemical:{}".format(len([item for values in self.dc.values() for item in values])))
        print("Disease-Mutation:{}".format(len([item for values in self.dm.values() for item in values])))
        print("Disease-Pathway:{}".format(len([item for values in self.dpa.values() for item in values])))
        print("Disease-Phenotype:{}".format(len([item for values in self.dpe.values() for item in values])))
        print("Disease-Disease:{}".format(len([item for values in self.dd.values() for item in values])))
        print("Gene-Chemical:{}".format(len([item for values in self.gc.values() for item in values])))
        print("Gene-Mutation:{}".format(len([item for values in self.gm.values() for item in values])))
        print("Gene-Pathway:{}".format(len([item for values in self.gpa.values() for item in values])))
        print("Gene-Phenotype:{}".format(len([item for values in self.gpe.values() for item in values])))
        print("Gene-Gene:{}".format(len([item for values in self.gg.values() for item in values])))
        print("Gene-Disease:{}".format(len([item for values in self.gd.values() for item in values])))
        print("Disease-Gene:{}".format(len([item for values in self.dg.values() for item in values])))

class BuildW():
    def __init__(self,chemical2id,mutation2id,pathway2id,phenotype2id,disease2indexpath,
                gene2indexpath,gene2id,disease2id,dc,dd,dm,dpa,dpe,gc,gg,gm,gpa,gpe,dg,gd):
        self.chemical2id=chemical2id
        self.mutation2id=mutation2id
        self.pathway2id=pathway2id
        self.phenotype2id=phenotype2id
        self.gene2id=gene2id
        self.disease2id=disease2id
        self.dc = dc
        self.dd=dd
        self.dm=dm
        self.dpa=dpa
        self.dpe=dpe
        self.gc=gc
        self.gg=gg
        self.gm=gm
        self.gpa=gpa
        self.gpe=gpe
        self.gd = gd
        self.dg = dg
        self.disease2indexpath = disease2indexpath
        self.gene2indexpath = gene2indexpath

        
    def buidw(self):
        """
        :return: csr格式
        """
        # for portrait we need rebuild data index
        # disease index+gene index
        maxdiseaseindex = max(self.disease2id.values())
        p_new_diseaseindex=self.disease2id
        p_new_geneindex = {gene:i+maxdiseaseindex+1 for gene,i in self.gene2id.items()}
        np.save(self.disease2indexpath,p_new_diseaseindex)
        np.save(self.gene2indexpath,p_new_geneindex)
        # for transverse we need rebuild data index
        # chemical index + disease index + mutation index + pathway index +phenotype index+ gene index
        t_new_chemicalindex = self.chemical2id
        maxchemicalindex=max(self.chemical2id.values())
        t_new_diseaseindex={disease:i+maxchemicalindex+1 for disease,i in self.disease2id.items()}# gene-disease-CNV ; disease-disease
        maxdiseaseindex = max(t_new_diseaseindex.values())
        t_new_mutationndex = {mutation:i+maxdiseaseindex+1 for mutation,i in self.mutation2id.items()}
        maxpathwayindex = max(t_new_mutationndex.values())
        t_new_pathwayindex = {pathway:i+maxpathwayindex+1 for pathway,i in self.pathway2id.items()}
        maxphenotype = max(t_new_pathwayindex.values())
        t_new_phenotypeindex = {phe:i+maxphenotype+1 for phe,i in self.phenotype2id.items()}
        maxgene = max(t_new_phenotypeindex.values())
        t_new_geneindex = {g:i+maxgene+1 for g,i in self.gene2id.items()}
        #构建CRS矩阵
        row=[]
        col=[]
        data=[]
        # build self.dc
        row, col, data = self.convert2csr(self.dc,p_new_diseaseindex,t_new_chemicalindex,row,col,data)
        # build self.dd
        row, col, data = self.convert2csr(self.dd, p_new_diseaseindex, t_new_diseaseindex, row, col, data)
        #build self.dm
        row, col, data = self.convert2csr(self.dm, p_new_diseaseindex, t_new_mutationndex, row, col, data)
        #build self.dpa
        row, col, data = self.convert2csr(self.dpa, p_new_diseaseindex, t_new_pathwayindex, row, col, data)
        # build self.dpe
        row, col, data = self.convert2csr(self.dpe, p_new_diseaseindex, t_new_phenotypeindex, row, col, data)
        # build self.gc
        row, col, data = self.convert2csr(self.gc, p_new_geneindex, t_new_chemicalindex, row, col, data)
        #build self.gg
        row, col, data = self.convert2csr(self.gg, p_new_geneindex, t_new_geneindex, row, col, data)
        #build self.gm
        row, col, data = self.convert2csr(self.gm, p_new_geneindex, t_new_mutationndex, row, col, data)
        #build self.gpa
        row, col, data = self.convert2csr(self.gpa, p_new_geneindex, t_new_pathwayindex, row, col, data)
        #build self.gpe
        row, col, data = self.convert2csr(self.gpe, p_new_geneindex, t_new_phenotypeindex, row, col, data)
        
        # Here if we use copuy number variation i the the net, then the gene-disease-relation will reused it 
        #build self.gd
        row, col, data = self.convert2csr(self.gd, p_new_geneindex, t_new_diseaseindex, row, col, data)
        #build self.dg
        row, col, data = self.convert2csr(self.dg, t_new_diseaseindex, p_new_geneindex, row, col, data)

        print("The shape of W: {}*{}".format(max(p_new_geneindex.values())+1,max(t_new_geneindex.values())+1))
        return csr_matrix((data,(row,col)),shape=(max(p_new_geneindex.values())+1,max(t_new_geneindex.values())+1))

    def convert2csr(self,relation,pindex,tindex,row,col,data,direction=False):
        for p,t in relation.items():
            for ti in t:
                row.append(pindex[p])
                if direction:
                    col.append(tindex[ti[0]])
                    data.append(int(ti[1]))
                    # data.append(1)
                else:
                    col.append(tindex[ti])
                    data.append(1)
        return row,col,data

class BuildX():
    def __init__(self,triplefile,disease2id,gene2id,outputX):
        """
        :param triplefile:
        :param relationship: how many relations we will use
        """
        self.path = triplefile
        self.disease2id=disease2id
        self.gene2id = gene2id
        self.outputX = outputX
        # self.function2id = function2id# LOF:0,GOF:1,COM:2
        maxdiseaseindex = max(self.disease2id.values())
        self.p_new_diseaseindex = self.disease2id
        self.p_new_geneindex = {gene: i + maxdiseaseindex + 1 for gene, i in self.gene2id.items()}
        self.maxp_new_geneindex = max(self.p_new_geneindex.values())

    def buildtensor(self):
        #here we rebuild index, these entity should has the same index as in the BuildW function.
        lofslice = self.get_slice(functionname="LOF")
        gofslice = self.get_slice(functionname="GOF")
        X = [lofslice,gofslice]
        print("The shape of X: {}*{}*{}".format(self.maxp_new_geneindex+1,self.maxp_new_geneindex+1,len(X)))
        np.save(self.outputX,X)
        return X

    def get_slice(self,functionname):
        rf = open(self.path,'r')
        row=[]
        col=[]
        data=[]
        pool = set()
        for line in rf:
            contents = line.strip().split()
            strfunction = contents[1]
            genes = contents[0].split(";")
            count = contents[3]
            for g in genes:
                gene = self.p_new_geneindex.get(g,-1)
                if gene!=-1:
                    disease = self.p_new_diseaseindex["MESH:"+contents[2]]
                    if (gene,disease) not in pool:
                        if strfunction==functionname or strfunction=="COM":
                            row.append(gene)
                            col.append(disease)
                            data.append(1)
                            # data.append(int(count))
                            row.append(disease)
                            col.append(gene)
                            data.append(1)
                            # data.append(int(count))
                            pool.add((gene,disease))
        return csr_matrix((data,(row,col)),shape=(max(self.p_new_geneindex.values())+1,max(self.p_new_geneindex.values())+1))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambdaA",help="parameter of lambdaA",type=float,required=True)
    parser.add_argument("--lambdaR",help="parameter of lambdaR",type=float,required=True)
    parser.add_argument("--lambdaV",help="parameter of lambdaV",type=float,required=True)
    parser.add_argument("--lambdaS",help="parameter of lambdaS",type=float,required=False)
    parser.add_argument("--maxiter",help="max iteration",type=int,default=100,required=False)
    parser.add_argument("--init_method",help="init method",type=str,default="nndsvd")
    parser.add_argument("--fit_method",help="fit_method",default=False)
    parser.add_argument("--checkpoint",help="checkpoint",type=str)

    parser.add_argument("--outputA",help="outputA",default="output/A.npy",type=str)
    parser.add_argument("--outputR",help="outputR",default="output/R.npy",type=str)
    parser.add_argument("--inputall",help="inputall",default=None,type=str)
    parser.add_argument("--inputtrain",help="inputtrain",default=None,type=str)
    parser.add_argument("--outputX",help="outputX",default="output/X.npy",type=str)
    parser.add_argument("--disease2indexpath",help="disease2indexpath",default="output/disease2index.npy",type=str)
    parser.add_argument("--gene2indexpath",help="gene2indexpath",default="output/gene2index.npy",type=str)
    parser.add_argument("--basedictionary",help="basedictionary",default="/mnt/disk2/kyzhou/毕业设计第二部分/DataForWA",type=str)
    parser.add_argument("--model",help="Use RESCAL or RESCAlSide",default="RESCALSide",type=str)
    
    args = parser.parse_args()
    maxiter = args.maxiter
    lambdaA = args.lambdaA
    lambdaR = args.lambdaR
    lambdaV = args.lambdaV
    lambdaS = args.lambdaS
    model = args.model
    fit_method = args.fit_method
    init_method = args.init_method
    checkpoint = args.checkpoint
    outputA = args.outputA
    outputR = args.outputR
    inputall = args.inputall
    inputtrain = args.inputtrain
    outputX = args.outputX
    disease2indexpath = args.disease2indexpath
    gene2indexpath = args.gene2indexpath
    basedictionary = args.basedictionary
    checkpoint = [int(x) for x in checkpoint.split(',')]

    # 将side information map到index上面去
    SM = SideMapping(basedictionary,inputall)
    chemical2id,mutation2id,pathway2id,phenotype2id,gene2id,disease2id = SM.entity2id()
    print("Chemical:{}\tMutation:{}\tPathway:{}\tPhenotype:{}"
          "\tGenes:{}\tDiseases:{}".format(len(chemical2id),len(mutation2id),len(pathway2id),len(phenotype2id),len(gene2id),len(disease2id)))
    SM.counter()# count some statistics
    #构建边缘矩阵W
    BW = BuildW(chemical2id,mutation2id,pathway2id,phenotype2id,disease2indexpath,
                gene2indexpath,gene2id,disease2id,SM.dc,SM.dd,SM.dm,SM.dpa,SM.dpe,
                 SM.gc,SM.gg,SM.gm,SM.gpa,SM.gpe,SM.dg,SM.gd)
    W = BW.buidw()# build matrix W
    #Build tensor X
    BA = BuildX(inputtrain,disease2id,gene2id,outputX)
    X = BA.buildtensor()
    #update
    if model=="RESCALSide":
        RE = RescalSide.RescalSide(maxiter,init_method,lambdaA,lambdaR,lambdaV,lambdaS,fit_method,checkpoint=checkpoint)
        A, R, V, iter, exctime = RE.rescalside(X,W,rank=100)
    elif model=="RESCAL":
        RE = Rescal.Rescal(maxiter,init_method,lambdaA,lambdaR,lambdaV,fit_method)
        A, R, f, iter, exctime = RE.rescal(X,rank=100)
    elif model=="PCA":
        RE = PCA.PCAN()
        A,R = RE.pca(W,dim=100)
    # the output of A will be a shape of 
    np.save(outputA,A)
    np.save(outputR,R)
if __name__=="__main__":
    main()
