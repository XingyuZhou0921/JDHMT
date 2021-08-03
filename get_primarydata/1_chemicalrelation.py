#! usr/bin/env python3
# -*- coding:utf-8 -*-
"""
@Author:zhoukaiyin
该脚本主要用于获取并处理
disease-drug
gene-drug
关系
主要数据来源是CTD
"""
import re
import numpy as np
from collections import defaultdict
def get_chemical_disease(filename,outputname):
    """
    chemical id use MESH
    disease id use MESH OR OMIM
    return chemical_id---disease_id
    """
    rf = open(filename,'r')
    wf = open(outputname,'w')
    wf.write("ChemicalID\tDiseaseID\tOmimIDs\tPubMedIDs\n")
    for line in rf:
        line = re.sub('"(.*?)"',"REPLACE",line)
        if line.startswith("#"):
            pass
        else:
            contents = line.strip().split(',')
            chemical_id = contents[1]
            disease_id = contents[4]
            OmimIDs = contents[8]
            PubMedIDs = contents[9]
            line = "{}\t{}\t{}\t{}\n".format(chemical_id,disease_id,OmimIDs,PubMedIDs)
            wf.write(line)
    rf.close()
    wf.close()

def chemical_gene(filename,outputname):
    rf = open(filename,'r')
    wf = open(outputname,'w')
    wf.write("ChemicalID\tGeneID\tOrganismID\tInteractionAction\tPubMedIDs\n")
    for line in rf:
        line = re.sub('"(.*?)"',"REPLACE",line)
        if line.startswith("#"):
            pass
        else:
            contents = line.strip().split(',')
            chemical_id = contents[1]
            gene_id = contents[4]
            InteractionActions = contents[9]
            PubMedIDs = contents[10]
            OrganismID = contents[7]
            line = "{}\t{}\t{}\t{}\t{}\n".format(chemical_id,gene_id,OrganismID,InteractionActions,PubMedIDs)
            wf.write(line)
    rf.close()
    wf.close()

def chemical_pathway(filename,outputname):
    def _norm(value,maxi,mini):
        return (value-mini)/(maxi-mini)
    wf = open(outputname,'w')
    wf.write("ChenmicalID\tPathwayID\tPvalue\n")
    chemical_pathway_pvalue = defaultdict(list)
    with open(filename) as rf:
        for line in rf:
            if line.startswith("#"):
                pass
            else:
                line = re.sub('"(.*?)"',"REPLACE",line)
                contents = line.split(',')
                chemicalid = contents[1]
                pathwayid = contents[4]
                try:
                    pvalue = -np.log(float(contents[5])+1e-100)
                except RuntimeWarning:
                    print(float(contents[5]))
                
                chemical_pathway_pvalue[chemicalid].append({pathwayid:pvalue})
    normal_chemical_pathway_pvalue = defaultdict(list)
    pvalues=[list(pp.values())[0] for chid,p_p in chemical_pathway_pvalue.items() for pp in p_p ] 
    max_ = np.max(pvalues)
    min_ = np.min(pvalues)
    normpvalue = [_norm(item,max_,min_) for item in pvalues]
    for chid,p_p in chemical_pathway_pvalue.items():
        for pp in p_p:

            pi = list(pp.values())[0]
            pathway = list(pp.keys())[0]
            wf.write("{}\t{}\t{}\n".format(chid,pathway,_norm(pi,max_,min_)))
    wf.close()


def rebuilt(lis):
    for i,item in enumerate(lis):
        kind = []
        if item=="[" or item=="]":
            kind.append(item)
    all_sentence = []
    start = []
    end = []
    for i,item in lis:
        if item=="[":
            start.append(i)
        elif item=="]":
            end.append(i)
    if kind==["[,],[,]"]:
        for s,e in zip(start,end):
            sentence = ' '.join(lis[s:e+1])
            all_sentence.append(sentence)
    elif kind=="[,[,],]":
        for s,e in zip(start,reversed(end)):
            sentence = ' '.join(lis[s:e+1])
            all_sentence.append(sentence)
    return all_sentence
        
        

def chemical_phenotype(filename,outputname):
    wf = open(outputname,'w')
    wf.write("\tPathwayID\tPvalue\n")
    ress = set()
    with open(filename) as rf:
        for line in rf:
            if line.startswith("#"):
                pass
            else:
                contents = line.split(',')
                organism = contents[7]
                if organism=="9606":
                    chemicalid = contents[1]
                    phenotype_name = contents[3]
                    phenotype_sentence = contents[8]
                    pathwayid = contents[4]
                    direction = contents[9]
                    direction = [item.split('^')[0] for item in direction.split("|") if item.split('^')[-1]=="phenotype"]
                    if len(direction)==1 and direction[0]=="phenotype":
                        pass
                    elif len(direction)==2:
                        
                        nphenotype_name = re.sub("\[","#[#",phenotype_name)
                        nphenotype_name = re.sub("\]","#]#",phenotype_name)
                        if len([nphenotype_name])!=1:
                            print(phenotype_name)
                            print(nphenotype_name)
                            sentences = rebuilt(nphenotype_name)
                            print(sentences)
                        

                    


def main():
    filename = "./chemical/CTD_chemicals_diseases.csv"#12 04 2019
    outputname = "./chemical/chemical_disease.txt"
    filename1 = "./chemical/CTD_chem_gene_ixns.csv"#12 04 2019 
    outputname1 = "./chemical/chemical_gene.txt"
    filename2 = "./chemical/CTD_chem_pathways_enriched.csv"
    outputname2 = "./chemical/chemical_pathway.txt"
    filename3 = "./chemical/CTD_chem_pheno.csv"
    outputname3 = "./chemical/chemical_phenotype.txt"

    # get_chemical_disease(filename,outputname)
    # chemical_gene(filename1,outputname1)
    # chemical_pathway(filename2,outputname2)
    chemical_phenotype(filename3,outputname3)

if __name__=="__main__":
    main()






