# build graph input with edgelist format
# edgelist format:node1 node2 <weight_float, optional>
# v1 版本将删除CNV 和microrna的使用
import json
import os
from collections import defaultdict
import itertools
class NodeMapping2id():
    def __init__(self,dictionary,alltriples,train_triples,edgoutput,labeloutput,filter_flag=True,unsigned=True):
        """
        disease_disease          gene_gene
        disease_chemical*        gene_chemical
        disease_mutation         gene_mutation
        disease_pathway          gene_pathway
        disease_phenotype        gene_phenotype
                    
                    disease_direction_microrna   -- not use
                    disease_CNV_gene             -- not use
                    *gene-disease

        *chemical-GO
        *chemical-pathway
        *chemical-phenotype

        nodes:
        disease,gene,chemical,mutation,pathway,phenotype,microrna
        TODO:Pathway-Pathway

        unsigned:边之间是否有方向属性，OpenNE中的方法不可以计算属性边
        """
        self.dictionary = dictionary
        self.alltriples = alltriples
        self.train_triples = train_triples
        self.edgoutputpath = os.path.join(dictionary,edgoutput)
        self.labeloutput = os.path.join(dictionary,labeloutput)
        self.dc = json.load(open(os.path.join(self.dictionary, "disease_chemical.json")))
        self.dd = json.load(open(os.path.join(self.dictionary, "disease_disease.json")))
        self.dm = json.load(open(os.path.join(self.dictionary, "disease_mutation.json")))
        self.dpa = json.load(open(os.path.join(self.dictionary, "disease_pathway.json")))
        self.dpe = json.load(open(os.path.join(self.dictionary, "disease_phe.json")))
        self.gg = json.load(open(os.path.join(self.dictionary, "gene_gene.json")))
        self.gm = json.load(open(os.path.join(self.dictionary, "gene_mutation.json")))
        self.gpa = json.load(open(os.path.join(self.dictionary, "gene_pathway.json")))
        self.gpe = json.load(open(os.path.join(self.dictionary, "gene_phe.json")))
        self.gc = json.load(open(os.path.join(self.dictionary, "gene_chemical.json")))

        # filter 
        # 我们的过滤策略，以目标疾病，和基因为中心，挑出直接相关的基因和疾病，
        # 再根据这些基因和疾病挑出相关的其他关系网络
        importantdiseases=importantgenes=[]
        if filter_flag:
            importantdiseases,importantgenes = self.get_important_gene_and_disease(self.alltriples,self.gg,self.dd)
            
        self.pipline(importantdiseases,importantgenes,filter_flag)

        data = self.dc.items()|self.dd.items()|self.dm.items()|self.dpa.items() |\
            self.dpe.items()|self.gc.items()|self.gg.items()|self.gm.items()|self.gpa.items()|self.gpe.items()

        self.node2id,self.id2node= self.mapallnode2id(data)
        self.convert2edgistlist(self.node2id,self.id2node,self.edgoutputpath,data,unsigned)
        self.get_node_labels(self.node2id,self.labeloutput )

    def get_important_gene_and_disease(self,triples,gg,dd):
        importantdiseases=set()
        importantgenes = set()
        rf = open(self.alltriples)
        for line in rf:
            contents = line.strip().split() 
            disease = contents[2]
            genes = contents[0].split(";")
            for gene in genes:
                importantgenes.add(gene)
                if gene in gg:
                    importantgenes = importantgenes | set(gg[gene])
            importantdiseases.add(disease)
            ndisease = "MESH:"+disease
            if ndisease in dd:
                importantdiseases = importantdiseases | set([d.split(":")[-1] for d in dd[ndisease]])        
        return importantdiseases,importantgenes

    def pipline(self,importantdiseases,importantgenes,filter_flag):
        if filter_flag:
            self.dc = {key:tuple(value) for key,value in self.dc.items() if key.split(":")[-1] in importantdiseases}
            self.dd = {key:tuple(value) for key,value in self.dd.items() if key.split(":")[-1] in importantdiseases}
            self.dm = {key:tuple(value) for key,value in self.dm.items() if key.split(":")[-1] in importantdiseases}
            self.dpa = {key:tuple(value) for key,value in self.dpa.items() if key.split(":")[-1] in importantdiseases}
            self.dpe = {key:tuple(value) for key,value in self.dpe.items() if key.split(":")[-1] in importantdiseases}
            self.gc = {key:tuple(value) for key,value in self.gc.items() if key in importantgenes}
            self.gg = {key:tuple(value) for key,value in self.gg.items() if key in importantgenes}
            self.gm = {key:tuple(value) for key,value in self.gm.items() if key in importantgenes}
            self.gpa = {key:tuple(value) for key,value in self.gpa.items() if key in importantgenes}
            self.gpe = {key:tuple(value) for key,value in self.gpe.items() if key in importantgenes}
            
        else:
            self.dc = {key:tuple(value) for key,value in self.dc.items()}
            self.dd = {key:tuple(value) for key,value in self.dd.items()}
            self.dm = {key:tuple(value) for key,value in self.dm.items()}
            self.dpa = {key:tuple(value) for key,value in self.dpa.items()}
            self.dpe = {key:tuple(value) for key,value in self.dpe.items()}
            self.gc = {key:tuple(value) for key,value in self.gc.items()}
            self.gg = {key:tuple(value) for key,value in self.gg.items()}
            self.gm = {key:tuple(value) for key,value in self.gm.items()}
            self.gpa = {key:tuple(value) for key,value in self.gpa.items()}
            self.gpe = {key:tuple(value) for key,value in self.gpe.items()}
            

    def get_node_labels(self,node2id,labeloutput):
        id2label = {}
        for node,id_ in node2id.items():
            if node.startswith("MESH:"):#disease
                id2label[id_]=1
            elif node.startswith("rs"):#mutation
                id2label[id_]=2
            elif node.startswith("HP"):# phenotype
                id2label[id_]=3
            elif node.isdigit():# gene
                id2label[id_]=4
            elif node.startswith("REACT") or node.startswith("KEGG"):# pathway
                id2label[id_]=5
            else:
                id2label[id_]=6 #chemical
        with open(labeloutput,'w') as wf:
            for id_,label in id2label.items():
                wf.write("{}\t{}\n".format(id_,label))

    def mapallnode2id(self,data):
        "nodes:disease,gene,chemical,mutation,pathway,phenotype"
        allnodes = set()
        "chenmical+'CHEM',mutation+'SNP:'"
        for key,value in data:
            allnodes.add(key)
            if isinstance(value[0],tuple):
                for item in value:
                    allnodes.add(item[0])
            else:
                allnodes.update(set(value))
        allnodes = sorted(list(allnodes))
        node2id= {item:i for i,item in enumerate(allnodes,1)}
        id2node= {i:item for i,item in enumerate(allnodes,1)}
        return node2id,id2node
    def train_data_edgist(self):
        train = set()
        with open(self.alltriples) as rf:
            for line in rf:
                contents = line.strip().split('\t')
                gene = contents[0].split(';')
                disease = contents[2]
                count = contents[3]
                for g in gene:
                    train.add((g,disease,count))
        return train


    def convert2edgistlist(self,node2id:dict,id2node:dict,output,data,unsigned=True):
        wf = open(output,'w')
        item2item = defaultdict(list)
        for item, listitem in data:
            for i in listitem:
                if isinstance(i,tuple):
                    node1 = node2id[item]
                    node2 = node2id[i[0]]
                    score = i[1]
                    if unsigned:
                         wf.write("{}\t{}\t{}\n".format(node1,node2,1))
                    else:
                        wf.write("{}\t{}\t{}\n".format(node1,node2,score))
                else:
                    node1 = node2id[item]
                    node2 = node2id[i]
                   
                    wf.write("{}\t{}\t{}\n".format(node1,node2,1))
        for g,d,c in self.train_data_edgist():
            try:
                # wf.write("{}\t{}\t{}\n".format(node2id[g],node2id["MESH:"+d],c))
                wf.write("{}\t{}\t{}\n".format(node2id[g],node2id["MESH:"+d],1))
            except KeyError:
                print(g)
            # return item2item

    def mapindex2nodes(self,outputs):
        with open(outputs,'w') as wf:
            wf.write("Index\tNodes\n")
            for index,node in self.id2node.items():
                wf.write("{}\t{}\n".format(index,node))
      
dictionary="./DataForWA"
alltriples="./triples/evidence.txt"
train_triples="/triples/evidence_train.txt"
edgoutput="side_ediges.txt"
labeloutput="nodelabels.txt"
filter_flag=False
# unsigned=True
id2nodes="DataForWA/id2nodes.txt"
Nodemap = NodeMapping2id(dictionary,alltriples,train_triples,edgoutput,labeloutput,filter_flag=False)
Nodemap.mapindex2nodes(id2nodes)
    


