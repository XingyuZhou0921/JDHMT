#ďź-*-encoding:utf-8-*-
import json
from collections import defaultdict
import numpy as np
from mapping import *
def get_diseases(file):
    diseases = set()
    with open(file) as rf:
        for line in rf:
            contents = line.strip().split('\t')
            diseases.add("MESH:"+contents[-1])
    return diseases

def get_genes(file):
    genes = set()
    with open(file) as rf:
        for line in rf:
            contents = line.strip().split('\t')
            genes.add(contents[-1])
    return genes

def build_disease_chemical(infile,outfile,diseases):
    """
    infile:chemical_disease.txt
    """
    disease_chemical = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.split('\t')
            chemical = contents[0]
            disease = contents[1]
            if disease in diseases:
                disease_chemical[disease].add(chemical)
    disease_chemical = {disease:list(chemical) for disease,chemical in disease_chemical.items()}    
    json.dump(disease_chemical,open(outfile,'w'))
    all_chemical = []
    all_disease = []
    for d,c in disease_chemical.items():
        all_disease.append(d)
        all_chemical.extend(c)
    return all_chemical,all_disease
    
def build_gene_chemical(infile,outfile,genes):
    gene_chemical = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.split('\t')
            chemical = contents[0]
            gene = contents[1]
            org = contents[2]
            if gene in genes and org in ["9606",""]:
                gene_chemical[gene].add(chemical)
    gene_chemical = {gene:list(chemical) for gene,chemical in gene_chemical.items()}   
    json.dump(gene_chemical,open(outfile,'w'))
    all_chemical = []
    all_gene = []
    for d,c in gene_chemical.items():
        all_gene.append(d)
        all_chemical.extend(c)
    return all_chemical,all_gene
   
        
def build_disease_pathway(infile,outfile,diseases):
    disease_pathway = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.strip().split('\t')
            pathway = contents[0]
            disease = contents[1]
            if disease in diseases:
                disease_pathway[disease].add(pathway)
    disease_pathway = {disease:list(pathway) for disease,pathway in disease_pathway.items()}   
    json.dump(disease_pathway,open(outfile,'w'))
    all_pathway = []
    all_disease = []
    for d,c in disease_pathway.items():
        all_disease.append(d)
        all_pathway.extend(c)
    return all_pathway,all_disease
    

def build_gene_pathway(infile,outfile,genes):
    gene_pathway = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.strip().split('\t')
            pathway = contents[0]
            gene = contents[1]
            if gene in genes:
                gene_pathway[gene].add(pathway)
    gene_pathway = {gene:list(pathway) for gene,pathway in gene_pathway.items()} 
    json.dump(gene_pathway,open(outfile,'w'))
    all_pathway = []
    all_gene = []
    for d,c in gene_pathway.items():
        all_gene.append(d)
        all_pathway.extend(c)
    return all_pathway,all_gene


def build_disease_phenotype(infile,outfile,diseases):
    disease_phenotype = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.strip().split('\t')
            pho = contents[0]
            disease = contents[1]
            if disease in diseases:
                disease_phenotype[disease].add(pho)
    disease_phenotype = {disease:list(phenotype) for disease,phenotype in disease_phenotype.items()} 
    json.dump(disease_phenotype,open(outfile,'w')) 
    all_phenotype = []
    all_disease = []
    for d,c in disease_phenotype.items():
        all_disease.append(d)
        all_phenotype.extend(c)
    return all_phenotype,all_disease

def build_gene_phenotype(infile,outfile,genes):
    gene_phenotype = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.strip().split('\t')
            pho = contents[0]
            gene = contents[1]
            if gene in genes:
                gene_phenotype[gene].add(pho)
    gene_phenotype = {gene:list(phenotype) for gene,phenotype in gene_phenotype.items()} 
    json.dump(gene_phenotype,open(outfile,'w'))
    all_phenotype = []
    all_gene = []
    for d,c in gene_phenotype.items():
        all_gene.append(d)
        all_phenotype.extend(c)
    return all_phenotype,all_gene
    
def build_gene_mutation(infile,outfile,genes):
    gene_mutation = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.strip().split('\t')
            mutation = contents[0]
            gene = contents[1]
            if gene in genes:
                gene_mutation[gene].add(mutation)
    gene_mutation = {gene:list(mutation) for gene,mutation in gene_mutation.items()} 
    json.dump(gene_mutation,open(outfile,'w'))
    all_mutation= []
    all_gene = []
    for d,c in gene_mutation.items():
        all_gene.append(d)
        all_mutation.extend(c)
    return all_mutation,all_gene

def build_disease_mutation(infile,outfile,diseases):
    disease_mutation = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.strip().split('\t')
            mu = contents[0]
            disease = contents[1]
            if disease in diseases:
                disease_mutation[disease].add(mu)
    disease_mutation = {disease:list(mu) for disease,mu in disease_mutation.items()} 
    json.dump(disease_mutation,open(outfile,'w')) 
    all_mutation= []
    all_disease = []
    for d,c in disease_mutation.items():
        all_disease.append(d)
        all_mutation.extend(c)
    return all_mutation,all_disease

def build_gene_gene(infile,outfile,genes):
    gene_gene = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.strip().split('\t')
            geneA = contents[0]
            geneB = contents[1]
            if geneA in genes and geneB in genes:
                gene_gene[geneA].add(geneB)
    gene_gene = {geneA:list(geneB) for geneA,geneB in gene_gene.items()} 
    json.dump(gene_gene,open(outfile,'w'))
    all_genes=[]
    for ga,gb in gene_gene.items():
        all_genes.append(ga)
        all_genes.extend(gb)
    return all_genes

def build_disease_disease(infile,outfile,diseases):
    disease_disease = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.strip().split('\t')
            disease1 = contents[0]
            disease = contents[1]
            if disease in diseases and disease1 in diseases:
                disease_disease[disease].add(disease1)
    disease_disease = {disease:list(disease1) for disease,disease1 in disease_disease.items()} 
    json.dump(disease_disease,open(outfile,'w')) 
    all_diseases=[]
    for da,db in disease_disease.items():
        all_diseases.append(da)
        all_diseases.extend(db)
    return all_diseases

def build_gene_disease(infile,outfile1,outfile2,genes,diseases):
    gene_disease = defaultdict(set)
    disease_gene = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.strip().split('\t')
            gene = contents[0]
            disease = contents[1]
            if gene in genes and disease in diseases:
                gene_disease[gene].add(disease)
                disease_gene[disease].add(gene)
    gene_disease = {gene:list(disease) for gene,disease in gene_disease.items()} 
    disease_gene = {disease:list(gene) for disease,gene in disease_gene.items()} 
    json.dump(gene_disease,open(outfile1,'w')) 
    json.dump(disease_gene,open(outfile2,'w')) 
    all_diseases=[]
    all_genes = []
    for g,d in gene_disease.items():
        all_genes.append(g)
        all_diseases.extend(d)
    return all_diseases,all_genes

def get_gene_cnv_disease(infile,outfile,genes,diseases):
    disease_cnv_gene = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.split('\t')
            gid = contents[0]
            cnv = contents[1]
            disease = contents[2].strip()
            disease_cnv_gene[disease].add((gid,cnv))
    disease_cnv_gene = {disease:list(gi_cn) for disease,gi_cn in disease_cnv_gene.items()} 
    json.dump(disease_cnv_gene,open(outfile,'w'))
    all_gene= []
    all_disease = []
    for d,c in disease_cnv_gene.items():
        all_disease.append(d)
        for item in c:
            gene = item[0]
            all_gene.append(gene)
    return all_gene,all_disease

def get_disease_direction_micro(infile,outfile,diseases):
    disease_dire_micro = defaultdict(set)
    with open(infile) as rf:
        rf.readline()
        for line in rf:
            contents = line.split('\t')
            disease = contents[0]
            dire = contents[1]
            micr = contents[2].strip()
            disease_dire_micro[disease].add((micr,dire))
    disease_dire_micro = {disease:list(micr_dire) for disease,micr_dire in disease_dire_micro.items()} 
    json.dump(disease_dire_micro,open(outfile,'w'))
    all_micro= []
    all_disease = []
    for d,c in disease_dire_micro.items():
        all_disease.append(d)
        for item in c:
            gene = item[0]
            all_micro.append(gene)
    return all_micro,all_disease


def main():
    diseases = get_diseases("vocabs/diseasename2mesh.txt")
    genes = get_genes("vocabs/genename2geneid.txt")
    
    all_chemical1,all_disease1 = build_disease_chemical("chemical/chemical_disease.txt","DataForWA/disease_chemical.json",diseases)
    all_chemical2,all_gene1 = build_gene_chemical("chemical/chemical_gene.txt","DataForWA/gene_chemical.json",genes)

    all_pathway1,all_disease2 = build_disease_pathway("pathway/pathway_disease.txt","DataForWA/disease_pathway.json",diseases)
    all_pathway2,all_gene2 = build_gene_pathway("pathway/pathway_gene.txt","DataForWA/gene_pathway.json",genes)
    
    all_phenotype1,all_disease3 = build_disease_phenotype("phenotype/phenotype_disease.txt","DataForWA/disease_phe.json",diseases)
    all_phenotype2,all_gene3 = build_gene_phenotype("phenotype/phenotype_gene.txt","DataForWA/gene_phe.json",genes)

    all_mutation1,all_gene4 = build_gene_mutation("mutation/mutation_gene.txt","DataForWA/gene_mutation.json",genes)
    all_mutation2,all_disease4 = build_disease_mutation("mutation/mutation_disease.txt","DataForWA/disease_mutation.json",diseases)

    all_gene5 = build_gene_gene("Gene_Gene/gene_gene.txt","DataForWA/gene_gene.json",genes)

    all_disease5 = build_disease_disease("Disease_Disease/disease_disease.txt","DataForWA/disease_disease.json",diseases)

    all_gene6,all_disease6 = get_gene_cnv_disease("CNV/gene_cnv_disease.txt","DataForWA/gene_cnv_disease.json",genes,diseases)

    all_micro,all_disease6 = get_disease_direction_micro("microrna/mesh2mcro.txt","DataForWA/mesh2mcro.json",diseases)

    # all_disease7,all_gene7 = build_gene_disease("gene-disease/gene_disease.txt","DataForWA/gene_disease.json",genes,diseases)
    all_disease7,all_gene7 = build_gene_disease("gene-disease/gene_disease.txt","DataForWA/gene_disease.json","DataForWA/disease_gene.json",genes,diseases)

    all_chemical1.extend(all_chemical2)
    all_pathway1.extend(all_pathway2)
    all_phenotype1.extend(all_phenotype2)
    all_mutation1.extend(all_mutation2)
    all_gene1.extend(all_gene2)
    all_gene1.extend(all_gene3)
    all_gene1.extend(all_gene4)
    all_gene1.extend(all_gene5)
    all_gene1.extend(all_gene6)
    all_gene1.extend(all_gene7)
    all_disease1.extend(all_disease2)
    all_disease1.extend(all_disease3)
    all_disease1.extend(all_disease4)
    all_disease1.extend(all_disease5)
    all_disease1.extend(all_disease6)
    all_disease1.extend(all_disease7)
    
    print("Use chemical:{}\nUse pathway:{}\nUse phenotype:{}\nUse mutation:{}\nUse gene:{}\nUse Disease:{}\nUse Micro:{}".format(
        len(set(all_chemical1)),
        len(set(all_pathway1)),
        len(set(all_phenotype1)),
        len(set(all_mutation1)),
        len(set(all_gene1)),
        len(set(all_disease1)),
        len(set(all_micro))
        ))

if __name__=="__main__":
    main()