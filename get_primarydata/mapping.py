#! -*-encoding:utf-8-*-
import pickle
import os
import codecs
import re
def symbol2id():
    symbol2id = {}
    datapath = "vocabs/9606_gene.txt"
    with open(datapath,'r') as rf:
        for line in rf:
            contents = line.split('\t')
            gene_id = contents[1]
            symbol = contents[2]
            synomys = contents[4].split('|')
            if symbol not in symbol2id:
                symbol2id[symbol] = gene_id
            for sy in synomys:
                if sy not in symbol2id:
                    symbol2id[sy]=gene_id
    return symbol2id

def disease_mapping(file):
    "disease_mappings.tsv come from DisGeNET" 
    UMLS2MESH = {}
    with open(file) as rf:
        for line in rf:
            line = line.strip().split('|')
            umls = line[0]
            type_ = line[2]
            if type_=="MSH":
                UMLS2MESH[umls]=line[3]
    return UMLS2MESH

def disease_DO2mesh(file):
    "diseaseontolgy.obo"
    DO2mesh = {}
    a_list=[]
    with open(file) as rf:
        for line in rf:
            line = line.strip()
            if line=="[Term]":
                if len(a_list)==2:
                    DO2mesh[a_list[0]]=a_list[1]
                a_list=[]
            line_c = line.split(':')
            if line_c[0]=="id":
                a_list.append("DOID:"+line_c[-1].strip())
            elif line_c[0]=="xref" and line_c[1]==" MESH":
                a_list.append("MESH:"+line_c[-1].strip())
    return DO2mesh

def sampleID2NCI(ComsicSample):
    sampleid2nci={}
    with open(ComsicSample,'r') as rf:
        rf.readline()
        for line in rf:
            contents = line.split()
            sampleID = contents[0]
            NCI = re.findall("C\d+",line)
            if len(NCI)!=0:
                sampleid2nci[sampleID]=NCI[0]
    return sampleid2nci

def nci2mesh(diseaseontolgy,targetcni2mesh):
    "diseaseontolgy.obo here we need add target nci2mesh"
    nci2mesh = {}
    meshlist = []
    ncilist = []
    with open(diseaseontolgy) as rf:
        for line in rf:
            line = line.strip()
            if line=="[Term]":
                if len(meshlist)>=1 and len(ncilist)>=1:
                    for n in ncilist:
                        for m in meshlist:
                            nci2mesh[n]=m
                meshlist = []
                ncilist = []
            line_c = line.split(':')
            if line_c[0]=="xref"and line_c[1]==" NCI":
                ncilist.append(line_c[-1].strip())
            elif line_c[0]=="xref" and line_c[1]==" MESH":
                meshlist.append(line_c[-1].strip())
        for nci,mesh in targetcni2mesh.items():
            if nci not in nci2mesh:
                nci2mesh[nci]=mesh
    return nci2mesh

