from mapping import *
import xlrd
from collections import defaultdict
def get_diseases(filename,outputname,disease2mesh):
    rf = open(filename,'r')
    wf = open(outputname,'w')
    a_list=[]
    # wf.write("DiseaseName\tMESH\n")
    for line in rf:
        line = line.strip()
        if line=="[Term]":
            if len(a_list)==2:
                line = "{}\t{}\n".format(a_list[0].lower(),'\t'.join(a_list[1:]))
                # line = "{}\t{}\n".format(a_list[1])
                wf.write(line)  
            a_list = []
        line_c = line.split(':')
        if line_c[0]=="name":
            a_list.append(line_c[-1].strip())
        elif line_c[0]=="xref" and line_c[1]==" MESH":
            a_list.append(line_c[-1])
    for disease,mesh in disease2mesh.items():
        wf.write("{}\t{}\n".format(disease,list(mesh)[0]))
    rf.close()
    wf.close()

def get_genes(filename,outputname):
    #TODO 这里gene的选取需要注意判断，是选取所有基因还是选取包含突变的基因
    symbol2id = {}
    with open(filename,'r') as rf:
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
    with open(outputname,'w') as wf:
        for gene,id_ in symbol2id.items():
            wf.write("{}\t{}\n".format(gene,id_))

def get_disease_mesh_cni(target_disease):
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
    diseasen2mesh={}
    cni2mesh = {}
    for disease,values in hasdir.items():
        diseasen2mesh[disease]=values.keys()
    for mesh_cnis in hasdir.values():
        mesh = list(mesh_cnis.keys())[0]
        cnis = mesh_cnis[mesh]
        for item in cnis:
            if len(item)!=0:
                cni2mesh[item]=mesh
    return diseasen2mesh,cni2mesh      

def main():
    # disease contain two file, one come from diseaseontolgy, another come from our target disease
    filename = "vocabs/diseaseontolgy.obo"
    filename_add = ""
    outputname = "vocabs/diseasename2mesh.txt"
    filename1 = "vocabs/9606_gene.txt"
    outputname1 = "vocabs/genename2geneid.txt"
    target_disease = "./target_disease.xlsx"
    diseasen2mesh,cni2mesh = get_disease_mesh_cni(target_disease)
    get_diseases(filename,outputname,diseasen2mesh) #这里构建的疾病会有重复，并且同一个疾病可能有两个mesh号
    get_genes(filename1,outputname1)


if __name__=="__main__":
    main()
