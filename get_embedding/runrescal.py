#! usr/bin/env python3
# -*- coding:utf-8 -*-
"""
@Author:zhoukaiyin
"""
from RescalSideV1 import RescalSide
from scipy.sparse import csr_matrix
import json
import os
import argparse

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
    def __init__(self,dictionary):
        self.dictionary = dictionary
        self.dc = json.load(open(os.path.join(self.dictionary, "disease_chemical.json")))
        self.dd = json.load(open(os.path.join(self.dictionary, "disease_disease.json")))
        self.dm = json.load(open(os.path.join(self.dictionary, "disease_mutation.json")))
        self.dpa = json.load(open(os.path.join(self.dictionary, "disease_pathway.json")))
        self.dpe = json.load(open(os.path.join(self.dictionary, "disease_phe.json")))
        self.gc = json.load(open(os.path.join(self.dictionary, "gene_chemical.json")))
        self.gg = json.load(open(os.path.join(self.dictionary, "gene_gene.json")))
        self.gm = json.load(open(os.path.join(self.dictionary, "gene_mutation.json")))
        self.gpa = json.load(open(os.path.join(self.dictionary, "gene_pathway.json")))
        self.gpe = json.load(open(os.path.join(self.dictionary, "gene_phe.json")))


    def hebing(self,diseasdir, genedir):
        entity = []
        for item in diseasdir.values():
            entity.extend(item)
        for item1 in genedir.values():
            entity.extend(item1)

        return sorted(set(entity))

    def merge_all(self,dc,dm,dpa,dpe,dd):
        entity = [k for d in [dc,dm,dpa,dpe,dd] for k in d.keys()]
        for item2 in dd.values():
            entity.extend(item2)
        return sorted(set(entity))

    def entity2id(self):
        chemicals= self.hebing(self.dc,self.gc)
        mutations = self.hebing(self.dm,self.gm)
        pathways = self.hebing(self.dpa,self.gpa)
        phenotypes = self.hebing(self.dpe,self.gpe)
        diseases = self.merge_all(self.dc,self.dm,self.dpa,self.dpe,self.dd)
        genes = self.merge_all(self.gc,self.gm,self.gpa,self.gpe,self.gg)
        # map 2 id
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

class BuildW():
    def __init__(self,chemical2id,mutation2id,pathway2id,phenotype2id,
                 gene2id,disease2id,dc,dd,dm,dpa,dpe,gc,gg,gm,gpa,gpe):
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

    def buidw(self):
        """
        :return: csr格式
        """
        # for portrait we need rebuild data index
        # disease index+gene index
        maxdiseaseindex = max(self.disease2id.values())
        p_new_diseaseindex=self.disease2id
        p_new_geneindex = {gene:i+maxdiseaseindex+1 for gene,i in self.gene2id.items()}
        # for transverse we need rebuild data index
        # chemical index + disease index + mutation index + pathway index +phenotype index+ gene index
        t_new_chemicalindex = self.chemical2id
        maxchemicalindex=max(self.chemical2id.values())
        t_new_diseaseindex={disease:i+maxchemicalindex+1 for disease,i in self.disease2id.items()}#
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
        print("The shape of W: {}*{}".format(max(p_new_geneindex.values())+1,max(t_new_geneindex.values())+1))
        return csr_matrix((data,(row,col)),shape=(max(p_new_geneindex.values())+1,max(t_new_geneindex.values())+1))

    def convert2csr(self,relation,pindex,tindex,row,col,data):
        for p,t in relation.items():
            for ti in t:
                row.append(pindex[p])
                col.append(tindex[ti])
                data.append(1)
        return row,col,data

class BuildA():
    def __init__(self,triplefile,disease2id,gene2id):
        """
        :param triplefile:
        :param relationship: how many relations we will use
        """
        self.path = triplefile
        self.disease2id=disease2id
        self.gene2id = gene2id
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
            gene = self.p_new_geneindex.get(contents[0],-1)
            if gene!=-1:
                disease = self.p_new_diseaseindex["MESH:"+contents[2]]
                if (gene,disease) not in pool:
                    if strfunction==functionname or strfunction=="COM":
                        row.append(gene)
                        col.append(disease)
                        data.append(1)
                        row.append(disease)
                        col.append(gene)
                        data.append(1)
                        pool.add((gene,disease))
        return csr_matrix((data,(row,col)),shape=(max(self.p_new_geneindex.values())+1,max(self.p_new_geneindex.values())+1))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lambdaA",help="parameter of lambdaA",type=float,required=True)
    parser.add_argument("--lambdaR",help="parameter of lambdaR",type=float,required=True)
    parser.add_argument("--lambdaV",help="parameter of lambdaV",type=float,required=True)
    parser.add_argument("--lambdaS",help="parameter of lambdaS",type=float,default=1,required=False)
    parser.add_argument("--maxiter",help="max iteration",type=int,default=100,required=False)
    parser.add_argument("--init_method",help="init method",type=str,default="nndsvd")
    args = parser.parse_args()
    maxiter = args.maxiter
    lambdaA = args.lambdaA
    lambdaR = args.lambdaR
    lambdaV = args.lambdaV
    lambdaS = args.lambdaS
    init_method = args.init_method


    #///////////////////////////////////////////////////////////////
    #这部分想法错误，有了关系对后，gene，疾病的id不用单独映射，可以直接从关系对中获取
    # 将gene和疾病映射到对应的ID上
    # diseasemapping = Mapping("../vocabs/diseasename2mesh.txt")
    # genemapping = Mapping("../vocabs/genename2geneid.txt")
    # disease2id = diseasemapping.convert2id()
    # gene2id = genemapping.convert2id()
    # print(disease2id)
    #///////////////////////////////////////////////////////////////
    # 将side information map到index上面去
    SM = SideMapping("../DataForWA")
    chemical2id,mutation2id,pathway2id,phenotype2id,gene2id,disease2id = SM.entity2id()
    print("Chemical:{}\tMutation:{}\tPathway:{}\tPhenotype:{}"
          "\tGenes:{}\tDiseases:{}".format(len(chemical2id),len(mutation2id),len(pathway2id),len(phenotype2id),len(gene2id),len(disease2id)))
    SM.counter()# count some statistics
    #构建边缘矩阵W
    BW = BuildW(chemical2id,mutation2id,pathway2id,phenotype2id,
                 gene2id,disease2id,SM.dc,SM.dd,SM.dm,SM.dpa,SM.dpe,
                 SM.gc,SM.gg,SM.gm,SM.gpa,SM.gpe)

    W = BW.buidw()# build matrix W

    #构建张量A
    BA = BuildA("../DataForWA/gene_function_disease.txt",disease2id,gene2id)
    X = BA.buildtensor()
    #迭代更新

    RE = RescalSide.RescalSide(maxiter,init_method,lambdaA,
                    lambdaR,lambdaV,lambdaS,fit_method = None)
    A, R, V, iter, exctime = RE.rescalside(X,W,rank=2000)
    #


if __name__=="__main__":
    main()
