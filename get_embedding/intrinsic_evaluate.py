# encoding:utf-8
import numpy as np
from sklearn.manifold import TSNE
import seaborn as sns
def load_item2index(disease2index,gene2index):
    disease2index = np.load(open(disease2index,'rb'),allow_pickle=True)
    gene2index = np.load(open(gene2index,'rb'),allow_pickle=True)
    ndir = {**disease2index.item(),**gene2index.item()}
    return disease2index,gene2index,ndir

def load_index2embedding(A):
    index2embed = {}
    A = np.load(open(A,'rb'))
    for i,embed in enumerate(A):
        index2embed[i]=embed
    return index2embed,A

def get_itemset(embedding,item2index):
    item2embed={}
    for item,index in item2index.item().items():
        item2embed[item] = embedding[index].tolist()
    return item2embed

def load_index_mapping(file):
    "/mnt/disk2/kyzhou/毕业设计第二部分/DataForWA/id2nodes.txt"
    index2nodes = {}
    with open(file) as rf:
        for line in rf:
            contents = line.strip().split('\t')
            #print(contents[1])
            #0 and 1 diaohuan
            index2nodes[contents[1]]=contents[0]
    return index2nodes

def get_embedding(vectorpath,indexmapfile):
    #load index 
    index2nodes = load_index_mapping(indexmapfile)
    embeddings = {}
    with open(vectorpath,'r') as rf:
        firstline = rf.readline()
        for line in rf:
            contents = line.strip().split(' ')
            i = contents[0]
            emb = contents[1:]
            node = index2nodes[i]
            if node.startswith("MESH") or node.isdigit():
                embeddings[node]=emb# embeddings[item] = embedding
    dim = firstline.split(' ')[1]
    return embeddings,int(dim)

def conver_toarray(embeddings):
    "only for openNE model"
    embeddingarray=[]
    names = []
    for key,value in embeddings.items():
        embeddingarray.append(value)
        names.append(key)
    return names,embeddingarray

def extract(names,A_SNE):
    disease2embed = {}
    gene2embed = {}
    for name,embedding in zip(names,A_SNE):
        if name.startswith("MESH"):
            disease2embed[name]=embedding
        elif name.isdigit():
            gene2embed[name]=embedding
    return disease2embed,gene2embed 

def t_SNE(A):
    X_embed = TSNE(n_components=2).fit_transform(A)
    return X_embed.embedding_

def write_csv(outcsv,disease2embed,gene2embed,target_disease,usefilter):
    diseases = get_target_disease(target_disease)
    useddisease = set()
    with open(outcsv,'w') as wf:
        wf.write("label,name,c1,c2\n")
        for key,value in disease2embed.items():
            if usefilter:
                if key in diseases:
                    useddisease.add(key)
                    wf.write("{},{},{}\n".format(0,key,','.join([str(item) for item in value]))) 
            else:
                wf.write("{},{},{}\n".format(0,key,','.join([str(item) for item in value]))) 
        for key,value in gene2embed.items():
            wf.write("{},{},{}\n".format(1,key,','.join([str(item) for item in value]))) 
    print("一共使用了{}个疾病。".format(len(useddisease)))

def write_concat(output,embedding,labels,g_d):
    with open(output,'w') as wf:
        wf.write("label,name,c1,c2\n")
        for embed,label,name in zip(embedding,labels,g_d):
            wf.write("{},{},{}\n".format(label,name,",".join([str(item) for item in embed])))

def get_target_disease(target_disease):
    disease = set()
    with open(target_disease) as rf:
        for line in rf:
            contents = line.strip().split('\t')
            id=contents[-1]
            disease.add("MESH:"+id)
    return disease

def concat_gene_disease(glod_data,embeddingarray,items2index):
    "/mnt/disk2/kyzhou/毕业设计第二部分/triples/evidence.txt"
    # index2name = {value,key for key,value in item2index.items()}
    tuple_id = {}
    concat_embedding = []
    labels = []
    g_d = []
    with open(glod_data,'r') as rf:
        for line in rf:
            contents = line.strip().split('\t')
            geneid = contents[0]
            diseaseid = "MESH:"+contents[2]
            label = contents[1]
            geneindex = items2index.get(geneid,-1)
            diseaseindex = items2index.get(diseaseid,-1)
            if geneindex!=-1 and diseaseindex!=-1:
                gbed = embeddingarray[geneindex]
                dbed = embeddingarray[diseaseindex]
                concat_embedding.append(np.concatenate([gbed,dbed]))
                labels.append(label)
                g_d.append("{}_{}".format(geneid,diseaseid))
    return concat_embedding,labels,g_d

def concat_gene_disease(glod_data,embeddingarray,items2index,concat_method):
    "/mnt/disk2/kyzhou/2021.6.15/毕业设计第二部分/triples/evidence.txt"
    # index2name = {value,key for key,value in item2index.items()}
    tuple_id = {}
    concat_embedding = []
    labels = []
    g_d = []
    with open(glod_data,'r') as rf:
        for line in rf:
            contents = line.strip().split('\t')
            geneid = contents[0]
            diseaseid = "MESH:"+contents[2]
            label = contents[1]
            geneindex = items2index.get(geneid,-1)
            diseaseindex = items2index.get(diseaseid,-1)
            if geneindex!=-1 and diseaseindex!=-1:
                gbed = embeddingarray[geneindex]
                dbed = embeddingarray[diseaseindex]
                if concat_method=="reduce":
                    concat_embedding.append(np.array(gbed,dtype=float)-np.array(dbed,dtype=float))
                elif concat_method=="concat":
                    concat_embedding.append(np.concatenate([gbed,dbed]))
                labels.append(label)
                g_d.append("{}_{}".format(geneid,diseaseid))
    return concat_embedding,labels,g_d

def main():
    "output: {concat_method}_{mode}_{filter}"
    mode="RotatE" #rescal, rescalonly, deepwalk.....
    dataindex = "9"
    # concat_method ="concat" #single
    # intrinsic evaluation contain: single, reduce and concat format.
    concat_method ="reduce" #concat, single ,reduce 
    usefilter = False
    if usefilter:
        filter_ = "filter"
    else:
        filter_ = "unfilter"
    glod_data="/home/kyzhou/2021.6.15/毕业设计第二部分/triples/evidence_train.txt"
    if mode in {"Ro"}:
        disease2indexpath = "output/disease2index{}.npy".format(dataindex)
        gene2indexpath = "output/gene2index{}.npy".format(dataindex)
        Apath = "output/A{}.npy".format(dataindex)
        disease2index,gene2index,items2index = load_item2index(disease2indexpath,gene2indexpath)
        index2embed,A = load_index2embedding(Apath) #{1:...,2:...}
        if concat_method in {"concat","reduce"}:
            concat_embedding,labels,g_d=concat_gene_disease(glod_data ,A ,items2index,concat_method)
            A_SNE = TSNE(n_components=2,verbose=1,random_state=0).fit_transform(concat_embedding)
            write_concat("t_SNE/{}_ {}_{}.csv".format(concat_method,mode,filter_),A_SNE,labels,g_d)
        elif concat_method=="single":
            A_SNE = TSNE(n_components=2,verbose=1,random_state=0).fit_transform(A)
            disease2embed = get_itemset(A_SNE,disease2index)
            gene2embed = get_itemset(A_SNE,gene2index)
            write_csv("t_SNE/{}_{}_{}.csv".format(concat_method,mode,filter_),disease2embed,gene2embed,"targetdisease.txt",usefilter)
    elif mode in {"deepwalk","line","node2vec","grarep","graphfactorization","RotatE"}:
        indexmapfile = "/home/kyzhou/2021.6.15/毕业设计第二部分/DataForWA/id2nodes1.txt"
        vectorpath = "/home/kyzhou/2021.6.15/RESCALSIDE/GraphEmbedding/OpenNE/output/vector_{}_Unfilter".format(mode)
        embeddings,dim = get_embedding(vectorpath,indexmapfile)
        names,embeddingarray = conver_toarray(embeddings)
        items2index = {key:value for value,key in enumerate(names)}
        if concat_method in {"concat","reduce"}:
            concat_embedding,labels,g_d=concat_gene_disease(glod_data ,embeddingarray ,items2index,concat_method)
            A_SNE = TSNE(n_components=2,verbose=1,random_state=0).fit_transform(concat_embedding)
            print(1)#
            write_concat("t_SNE/{}_ {}_{}.csv".format(concat_method,mode,filter_),A_SNE,labels,g_d)
        elif concat_method=="single":
            A_SNE = TSNE(n_components=2,verbose=1,random_state=0).fit_transform(embeddingarray)
            disease2embed,gene2embed  = extract(names,A_SNE)
            write_csv("t_SNE/{}_{}_{}.csv".format(concat_method,mode,filter_),disease2embed,gene2embed,"targetdisease.txt",usefilter)

if __name__=="__main__":
    main()
