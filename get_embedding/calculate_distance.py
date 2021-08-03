# encoding:utf-8
# for clarify the results of our model here we calculate the distance between target mean cancer disease distance. with  
import numpy as np
from sklearn.manifold import TSNE
import seaborn as sns
from sklearn import metrics
from sklearn.metrics import pairwise_distances
import math

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
            index2nodes[contents[0]]=contents[1]
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
            disease2embed[name]=np.array(embedding,dtype='float')
        elif name.isdigit():
            gene2embed[name]=np.array(embedding,dtype='float')
    return disease2embed,gene2embed 

def t_SNE(A):
    X_embed = TSNE(n_components=2).fit_transform(A)
    return X_embed.embedding_

def single_distance(disease2embed,gene2embed,target_disease,usefilter):
# def write_csv(outcsv,disease2embed,gene2embed,target_disease,usefilter):
    diseases = get_target_disease(target_disease)
    # useddisease = set()
    # with open(outcsv,'w') as wf:
    #     wf.write("label,name,c1,c2\n")
    #     for key,value in disease2embed.items():
    #         if usefilter:
    #             if key in diseases:
    #                 useddisease.add(key)
    #                 wf.write("{},{},{}\n".format(0,key,','.join([str(item) for item in value]))) 
    #         else:
    #             wf.write("{},{},{}\n".format(0,key,','.join([str(item) for item in value]))) 
    #     for key,value in gene2embed.items():
    #         wf.write("{},{},{}\n".format(1,key,','.join([str(item) for item in value]))) 
    # print("一共使用了{}个疾病。".format(len(useddisease)))
    embedding = []
    labels = []
    for key,value in disease2embed.items():
        if usefilter:
            if key in diseases:
                labels.append(0)
                embedding.append(value)
        else:
            labels.append(0)
            embedding.append(value)

    for key,value in gene2embed.items():
        labels.append(1)
        embedding.append(value)
    re1 = metrics.silhouette_score(embedding, labels, metric='euclidean')
    #Peter J. Rousseeuw (1987). “Silhouettes: a Graphical Aid to the Interpretation and Validation of Cluster Analysis”. Computational and Applied Mathematics 20: 53–65. doi:10.1016/0377-0427(87)90125-7.
    re2 = metrics.calinski_harabasz_score(embedding, labels)
    # Caliński, T., & Harabasz, J. (1974). “A Dendrite Method for Cluster Analysis”. Communications in Statistics-theory and Methods 3: 1-27. doi:10.1080/03610927408827101.

    return re1,re2

def concat_distance(embedding,labels):
# def concat_Silhouette_Coefficient(output,embedding,labels,g_d):
    re1 = metrics.silhouette_score(embedding, labels, metric='euclidean')
    re2 = metrics.calinski_harabasz_score(embedding, labels)
    return re1,re2 
    # with open(output,'w') as wf:
    #     wf.write("label,name,c1,c2\n")
    #     for embed,label,name in zip(embedding,labels,g_d):
    #         wf.write("{},{},{}\n".format(label,name,",".join([str(item) for item in embed])))

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
                if concat_method=="reduce":
                    concat_embedding.append(np.array(gbed,dtype=float)-np.array(dbed,dtype=float))
                elif concat_method=="concat":
                    concat_embedding.append(np.concatenate([gbed,dbed]))
                labels.append(label)
                g_d.append("{}_{}".format(geneid,diseaseid))
    return concat_embedding,labels,g_d

def get_most_similarity(inputtarget,disease2embed,gene2embed,topK=10):
    disease_distance = {}
    gene_distance = {}
    embed1 = disease2embed[inputtarget]
    for disease, embed2 in disease2embed.items():
        if inputtarget!=disease:
            # simliarity = np.dot(embed1,embed2)/(np.linalg.norm(embed1)*(np.linalg.norm(embed2)))
            simliarity = np.sqrt(np.sum(np.square(np.array(embed1)-np.array(embed2))))
            if not math.isnan(simliarity):
                disease_distance[disease] = simliarity 
    red = sorted(disease_distance.items(),key=lambda x:x[1])
    result_disease = red[:topK]

    for gene,embed3 in gene2embed.items():
        # simliarity = np.dot(embed1,embed3)/(np.linalg.norm(embed1)*(np.linalg.norm(embed3)))
        simliarity = np.sqrt(np.sum(np.square(np.array(embed1)-np.array(embed3))))
        if not math.isnan(simliarity):
            gene_distance[gene] = simliarity 
    reg = sorted(gene_distance.items(),key=lambda x:x[1])
    result_gene = reg[:topK]
    return result_disease,result_gene

def mapmesh2diseasename(orgindata):
    mesh2diseasename = {}
    with open(orgindata) as rf:
        for line in rf:
            contents = line.strip().split("\t")
            mesh = "MESH:"+contents[1]
            name = contents[0]
            mesh2diseasename[mesh]=name
    return mesh2diseasename

def mapgeneid2name(orgindata):
    gid2name = {}
    with open(orgindata) as rf:
        for line in rf:
            contents = line.strip().split("\t")
            id_ = contents[1]
            name = contents[0]
            gid2name[id_]=name
    return gid2name

def write_format(inputdisease,sim_diseases,sim_genes,mesh2name,geneid2name):
    inputname = mesh2name[inputdisease]
    simdiseases = []
    simgenes = []
    for mesh,score in sim_diseases:
        name = mesh2name[mesh]
        simdiseases.append(name)
    for gene,score in sim_genes:
        gene = geneid2name.get(gene,"None")
        simgenes.append(gene)
    print("Similarity Disease:")
    print("{}:\n{}\n".format(inputname,'\n'.join(simdiseases)))
    print("Similarity Gene:")
    print("{}\n".format('\n'.join(simgenes)))
    



def main():
    "output: {concat_method}_{mode}_{filter}"
    li = [("rescal","11"),
          ("rescalonly","27"),
          ("SVD","9"),
          ("deepwalk","9"),
          ("line","9"),
          ("node2vec","9"),
          ("grarep","9"),
          ("graphfactorization","9")]
    inputdisease = "MESH:D014062"# #D001943 乳腺癌，D006394 阿兹海默症，D009404 肾病综合征,D007680,D001749 膀胱癌，D001932脑癌,D014594子宫瘤，D012004直肠癌，D008113肝癌,D015459白血病，癫痫D004827，舌癌D014062
    for mode,dataindex in li:
        print("{}\t{}:".format(mode,inputdisease))
        # mode="rescal" #rescal, rescalonly, deepwalk.....
        # dataindex = "11"
        
        """
        dataindex for different disease:
        rescal:11
        rescalonly:27
        SVD:9
        deepwalk:
        LINE:
        """
        # concat_method ="concat" #single
        # intrinsic evaluation contain: single, reduce and concat format.
        concat_method ="single" #concat, single ,reduce 
        usefilter = False
    
        glod_data="/mnt/disk2/kyzhou/毕业设计第二部分/triples/evidence.txt"
        diseas_mesh = "/mnt/disk1/RESCALSIDE/DataForWA/diseasename2mesh.txt"
        targetdisease = "/mnt/disk1/RESCALSIDE/RESCAL_SIDEINFORMATION/targetdisease.txt"
        gene2id = "/mnt/disk1/RESCALSIDE/DataForWA/genename2geneid.txt"
        meshdisease1 = mapmesh2diseasename(diseas_mesh)
        meshdisease2 = mapmesh2diseasename(targetdisease)
        meshdisease = {**meshdisease1,**meshdisease2}
        geneid2name = mapgeneid2name(gene2id)
        if mode in {"rescal","rescalonly","SVD"}:
            disease2indexpath = "output/disease2index{}.npy".format(dataindex)
            gene2indexpath = "output/gene2index{}.npy".format(dataindex)
            Apath = "output/A{}.npy".format(dataindex)
            disease2index,gene2index,items2index = load_item2index(disease2indexpath,gene2indexpath)
            index2embed,A = load_index2embedding(Apath) #{1:...,2:...}
            if concat_method in {"concat","reduce"}:
                concat_embedding,labels,g_d=concat_gene_disease(glod_data ,A ,items2index,concat_method)
                A_SNE = TSNE(n_components=2,verbose=1,random_state=0).fit_transform(concat_embedding)
                # write_concat("t_SNE/{}_ {}_{}.csv".format(concat_method,mode,filter_),A_SNE,labels,g_d)
                re1,re2 = concat_distance(concat_method,labels)
            elif concat_method=="single":
                # A_SNE = TSNE(n_components=2,verbose=1,random_state=0).fit_transform(A)
                # disease2embed = get_itemset(A_SNE,disease2index)
                # gene2embed = get_itemset(A_SNE,gene2index)
                # A_SNE = TSNE(n_components=2,verbose=1,random_state=0).fit_transform(A)
                disease2embed = get_itemset(A,disease2index)
                gene2embed = get_itemset(A,gene2index)
                # write_csv("t_SNE/{}_{}_{}.csv".format(concat_method,mode,filter_),disease2embed,gene2embed,"targetdisease.txt",usefilter)
                # re1,re2 = single_distance(disease2embed,gene2embed,"targetdisease.txt",usefilter)
                result_disease,result_gene = get_most_similarity(inputdisease,disease2embed,gene2embed,topK=10)
                write_format(inputdisease,result_disease,result_gene,meshdisease,geneid2name)
        elif mode in {"deepwalk","line","node2vec","grarep","graphfactorization"}:
            indexmapfile = "/mnt/disk2/kyzhou/毕业设计第二部分/DataForWA/id2nodes.txt"
            vectorpath = "/mnt/disk1/RESCALSIDE/GraphEmbedding/OpenNE/output/vector_{}_Unfilter".format(mode)
            embeddings,dim = get_embedding(vectorpath,indexmapfile)
            names,embeddingarray = conver_toarray(embeddings)
            items2index = {key:value for value,key in enumerate(names)}
            if concat_method in {"concat","reduce"}:
                concat_embedding,labels,g_d=concat_gene_disease(glod_data ,embeddingarray ,items2index,concat_method)
                A_SNE = TSNE(n_components=2,verbose=1,random_state=0).fit_transform(concat_embedding)
                # write_concat("t_SNE/{}_ {}_{}.csv".format(concat_method,mode,filter_),A_SNE,labels,g_d)
                re1,re2 = concat_distance(concat_method,labels)
            elif concat_method=="single":
                # A_SNE = TSNE(n_components=2,verbose=1,random_state=0).fit_transform(embeddingarray)
                # disease2embed,gene2embed  = extract(names,A_SNE)
                # A_SNE = TSNE(n_components=2,verbose=1,random_state=0).fit_transform(embeddingarray)
                disease2embed,gene2embed  = extract(names,embeddingarray)
                # write_csv("t_SNE/{}_{}_{}.csv".format(concat_method,mode,filter_),disease2embed,gene2embed,"targetdisease.txt",usefilter)
                # re1,re2 = single_distance(disease2embed,gene2embed,"targetdisease.txt",usefilter)
                result_disease,result_gene=get_most_similarity(inputdisease,disease2embed,gene2embed,topK=10)
                write_format(inputdisease,result_disease,result_gene,meshdisease,geneid2name)
        # print("Model:{}\nSilhouette Coefficient Results Is: {}\n Calinski-Harabasz Result Is: {}".format(mode,re1,re2))

if __name__=="__main__":
    main()
