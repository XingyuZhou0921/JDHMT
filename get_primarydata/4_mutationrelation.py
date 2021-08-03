#! usr/bin/env python3
# -*- coding:utf-8 -*-
"""
@Author:zhoukaiyin
该脚本主要用于获取并处理
disease-mutation
gene-mutation
"""
import re
from mapping import *

def get_mutation_phenotype(filename,outputname):
    """
    phenotype not equal disease
    """
    rf = open(filename,'r')
    wf = open(outputname,'w')
    rf.readline()
    wf.write("MutationID\tPhenotypeID\n")
    for line in rf:
        contents = line.strip().split('\t')
        mutation_id = contents[-1]
        phenotype_id = contents[12]
        line = "{}\t{}\n".format(mutation_id,phenotype_id)
        wf.write(line)
    rf.close()
    wf.close()

def get_mutation_gene(filename,outputname):
    rf = open(filename,'r')
    wf = open(outputname,'w')
    rf.readline()
    wf.write("MutationID\tGeneID\n")
    for line in rf:
        contents = line.strip().split('\t')
        mutation_id = contents[-1]
        gene_id = contents[3]
        line = "{}\t{}\n".format(mutation_id,gene_id)
        wf.write(line)
    rf.close()
    wf.close()

def disgene_mutation_disease(filename,outputname,umls2mesh):
    rf = open(filename,'r')
    wf = open(outputname,'w')
    rf.readline()
    wf.write("SNPID\tDiseaseID\n")
    for line in rf:
        contents = line.strip().split('\t')
        snpid = contents[0]
        umls = contents[5]
        distype=contents[7]
        try:
            line = "{}\tMESH:{}\n".format(snpid,umls2mesh[umls])
            wf.write(line)
        except KeyError:
            # print("{}\n".format(umls))
            pass
    rf.close()
    wf.close()

def disgen_mutation_gene(filename,outputname,umls2mesh):
    """
    """
    rf = open(filename,'r')
    wf = open(outputname,'w')
    rf.readline()
    wf.write("SNPID\tGeneID\n")
    for line in rf:
        contents = line.strip().split('\t')
        snpid = contents[0]
        gene_id = contents[1]
        line = "{}\t{}\n".format(snpid,gene_id)
        wf.write(line)
    rf.close()
    wf.close()


def main():
    filename = "./mutation/curated_variant_disease_associations.tsv"# database:DisGeNET
    #The file contains variant-disease associations from UNIPROT, ClinVar, GWASdb, and the GWAS catalog.
    outputname = "./mutation/mutation_disease.txt"

    filename1 = "./mutation/variant_to_gene_mappings.tsv"
    #The file contains the mappings of DisGeNET variants (dbSNP Identifiers) to NCBI Entrez identifiers according to dbSNP database
    outputname1 = "./mutation/mutation_gene.txt"

    filename2 = "./vocabs/disease_mappings.tsv" # DisGNET
    umls2mesh = disease_mapping(filename2)
    disgene_mutation_disease(filename,outputname,umls2mesh)
    disgen_mutation_gene(filename1,outputname1,umls2mesh)

if __name__=="__main__":
    main()






