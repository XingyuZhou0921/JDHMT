#! -*- encoding:utf-8 -*-
import re
import mapping
import os
import pickle
import xlrd
from collections import defaultdict
def get_disease_mesh_cni(target_disease):
    "here we need keep all the target disease in the dictionary"
    data = xlrd.open_workbook(target_disease)
    data.sheet_names()
    table = data.sheet_by_name('target_cancer')
    rows = table.nrows
    cols = table.ncols
    hasdir = defaultdict(dir)
    for rownum in range(1,rows):
        row_values = [item.strip() for item in table.row_values(rownum)]
        diseasename = row_values[0]
        mesh = row_values[1]
        cni = [item for item in row_values[2:] if item!="-" and item!='']
        if diseasename!="-" and mesh!="-":
            hasdir[diseasename]={mesh:cni}
    cni2mesh = {}
    for mesh_cnis in hasdir.values():
        mesh = list(mesh_cnis.keys())[0]
        cnis = mesh_cnis[mesh]
        for item in cnis:
            if len(item)!=0:
                cni2mesh[item]=mesh
    return cni2mesh

def gene2function2disease(cnvfile,gene_cnv_disease,nci2mesh,sampleid2nci,symbol2id):

    wf = open(gene_cnv_disease,'w')
    wf.write("GeneID\tCNV\tDiseaseID\n")
    gene2function2diseae = {}
    with open(cnvfile) as rf:
        rf.readline()
        for line in rf:
            contents = line.split()
            gene = contents[2].split('_')[0]
            cnv = contents[-4]
            sampleid = contents[3]
            # id_tumor.add(id_tumour)
            nci_id = sampleid2nci.get(sampleid,-1)
            gid = symbol2id.get(gene,-1)
            if gid!=-1 and nci_id!=-1 :
                mesh = nci2mesh.get(nci_id,-1)
                if mesh!=-1:
                    cnv = contents[-4]
                    if cnv=="gain":
                        cnv=-1
                    elif cnv=="loss":
                        cnv=1
                    gene2function2diseae[(gid,mesh)]=cnv
    for gid_mesh,fun in gene2function2diseae.items():
        wf.write("{}\t{}\t{}\n".format(gid_mesh[0],fun,"MESH:"+gid_mesh[1]))
    
def main():
    diseaseontloy = "vocabs/diseaseontolgy.obo"
    ComsicSample = "CNV/CosmicSample.tsv"
    cnvfile = "CNV/CosmicCompleteCNA.tsv"
    gene_cnv_disease = "CNV/gene_cnv_disease.txt"
    target_disease = "target_disease.xlsx"
    sampleid2nci = mapping.sampleID2NCI(ComsicSample)
    targetcni2mesh = get_disease_mesh_cni(target_disease)
    nci2mesh = mapping.nci2mesh(diseaseontloy,targetcni2mesh)
    symbol2id = mapping.symbol2id()
    gene2function2disease(cnvfile,gene_cnv_disease,nci2mesh,sampleid2nci,symbol2id)

if __name__=="__main__":
    main()
