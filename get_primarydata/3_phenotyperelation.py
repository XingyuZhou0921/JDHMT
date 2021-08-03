#! usr/bin/env python3
# -*- coding:utf-8 -*-
"""
@Author:zhoukaiyin
该脚本主要用于获取并处理
disease-phenotype
gene-phenotype
关系
主要数据来源是HPO
"""
import re
def get_phenotype_gene(filename,outputname):
    """
    gene id: entrez gene id
    phenotype:HPO
    """
    rf = open(filename,'r')
    wf = open(outputname,'w')
    wf.write("PhenotypeID\tGeneID\n")
    for line in rf:
        if line.startswith("#"):
            pass
        else:
            contents = line.strip().split('\t')
            geneid = contents[0]
            phrnotype_id = contents[-1]
            line = "{}\t{}\n".format(phrnotype_id,geneid)
            wf.write(line)
    rf.close()
    wf.close()

def get_phenotype_disease(filename,outputname):
    """
    """
    rf = open(filename,'r')
    wf = open(outputname,'w')
    a_list=[]
    wf.write("phenotypeID\tDiseaseID\n")  
    for line in rf:
        line = line.strip()
        if line=="[Term]":
            if len(a_list)==2:
                line = "{}\t{}\n".format(a_list[0],a_list[1])
                wf.write(line)  
            a_list = []
        line_c = line.split(':')
        if line_c[0]=="id":
            a_list.append("HP:"+line_c[-1])
        elif line_c[0]=="xref" and line_c[1]==" MSH":
            a_list.append("MESH:"+line_c[-1])
    rf.close()
    wf.close()

def main():
    filename = "./phenotype/genes_to_phenotype.txt"
    outputname = "./phenotype/phenotype_gene.txt"
    filename1 = "./phenotype/hp.obo"
    outputname1 = "./phenotype/phenotype_disease.txt"
    get_phenotype_gene(filename,outputname)
    get_phenotype_disease(filename1,outputname1)

if __name__=="__main__":
    main()






