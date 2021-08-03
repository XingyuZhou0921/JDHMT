#! usr/bin/env python3
# -*- coding:utf-8 -*-
"""
@Author:zhoukaiyin
该脚本主要用于获取并处理
disease-pathway
gene-pathway
关系
主要数据来源是CTD
"""
import re
def get_pathway_disease(filename,outputname):
    """
     id use MESH
    disease id use MESH OR OMIM
    return chemical_id---disease_id
    """
    rf = open(filename,'r')
    wf = open(outputname,'w')
    wf.write("PathwayID\tDiseaseID\n")
    for line in rf:
        line = re.sub('"(.*?)"',"REPLACE",line)
        if line.startswith("#"):
            pass
        else:
            contents = line.strip().split(',')
            disease_id = contents[1]
            pathway_id = contents[3]
            line = "{}\t{}\n".format(pathway_id,disease_id)
            wf.write(line)
    rf.close()
    wf.close()

def pathway_gene(filename,outputname):
    """
    GeneID NCBI gene identifier
    Pathway ID  KEGG or RECACTOME identifier
    """
    rf = open(filename,'r')
    wf = open(outputname,'w')
    wf.write("PathwayID\tGeneID\n")
    for line in rf:
        line = re.sub('"(.*?)"',"REPLACE",line)
        if line.startswith("#"):
            pass
        else:
            contents = line.strip().split(',')
            gene_id = contents[1]
            pathway_id = contents[3]
            line = "{}\t{}\n".format(pathway_id,gene_id)
            wf.write(line)
    rf.close()
    wf.close()

def main():
    filename = "./pathway/CTD_diseases_pathways.csv"
    outputname = "./pathway/pathway_disease.txt"
    filename1 = "./pathway/CTD_genes_pathways.csv"
    outputname1 = "./pathway/pathway_gene.txt"
    get_pathway_disease(filename,outputname)
    pathway_gene(filename1,outputname1)

if __name__=="__main__":
    main()






