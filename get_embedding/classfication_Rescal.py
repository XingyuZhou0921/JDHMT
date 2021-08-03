import numpy as np
from collections import defaultdict,Counter
from sklearn.ensemble import RandomForestClassifier,BaggingClassifier,AdaBoostClassifier,VotingClassifier
from sklearn import svm
from sklearn.metrics import precision_score,recall_score,f1_score
from sklearn.model_selection import cross_validate
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier
from tqdm import tqdm
import random
import argparse
import warnings

warnings.filterwarnings("ignore")
def labelmap():
    return {"LOF":0,"GOF":1,"COM":2}
def dumptrain_or_total(inputfile):
    label2index = labelmap()
    rf = open(inputfile,'r')
    hashdir = defaultdict(list)
    for line in rf:
        contents = line.strip().split('\t')
        genes = contents[0].split(';')
        for g in genes:
            func = contents[1]
            mesh = "MESH:"+contents[2]
            tuple_ = tuple(sorted([g,mesh]))
            hashdir[tuple_].append(func)
    for t,l in hashdir.items():
        if len(l)>=2:
            hashdir[t]=label2index["COM"]
        else:
            hashdir[t]=label2index[l[0]]
    return hashdir 

def get_embedding(A,disease2index,gene2index):
    #load index 
    p_new_diseaseindex = np.load(disease2index,allow_pickle=True)
    p_new_geneindex = np.load(gene2index,allow_pickle=True)
    index2disease = {value:key for key,value in p_new_diseaseindex.item().items()}
    index2gene = {value:key for key,value in p_new_geneindex.item().items()}
    index2item = {**index2disease,**index2gene}
    embeddings = {}
    for i,emb in enumerate(A):
        embeddings[index2item[i]]=emb
    return embeddings

def random_sample_negtive(totaldir,targetY):
    totalgene = {key[0] for key in totaldir.keys()}
    totaldisease = {key[1] for key in totaldir.keys()}
    # train negative
    negative = {}
    totalNegativefortrain = max(Counter(targetY).items(),key=lambda x:x[1])[1]
    num=0
    while num<totalNegativefortrain+20:
        gene = random.choice(list(totalgene))
        disease = random.choice(list(totaldisease))
        tup = (gene,disease)
        if tup not in totaldir:
            if tup not in negative:
                negative[tup]=len(set(targetY))
                num+=1
    return negative

def concateembedding(embedings,gene_meshdir,empty):
    X = []
    Y = []
    for g_mesh,label in gene_meshdir.items():
        gene = embedings.get(g_mesh[0],-1)
        mesh = embedings.get(g_mesh[1],-1)
        if not isinstance(gene,int) and not isinstance(mesh,int):
            X.append(np.concatenate((gene,mesh)))
            Y.append(label)
    return X,Y

def load_data(intputall,inputtrain,outputA,disease2index,gene2index):
    totaldir = dumptrain_or_total(intputall)
    traindir = dumptrain_or_total(inputtrain)
    # traindir = dumptrain_or_total("../DataForWA/gene_function_disease_total_done.txt")
    testdir = {key:value for key,value in totaldir.items()-traindir.items()}
    A = np.load(outputA,allow_pickle=True)
    embedings = get_embedding(A,disease2index,gene2index)
    empty = np.zeros_like(A[0])
    trainX,trainY = concateembedding(embedings,traindir,empty)
    testX,testY = concateembedding(embedings,testdir,empty)
    negative_train = random_sample_negtive(totaldir,trainY)
    negative_test = random_sample_negtive(totaldir,testY)
    n_trainX,n_trainY = concateembedding(embedings,negative_train,empty)
    n_testX,n_testY = concateembedding(embedings,negative_test,empty)
    trainX.extend(n_trainX)
    trainY.extend(n_trainY)
    testX.extend(n_testX)
    testY.extend(n_testY)
    trainX,trainY=_shuffle(trainX,trainY)
    testX,testY=_shuffle(testX,testY)
    return trainX,trainY,testX,testY
    
def _shuffle(trainX,trainY):
    n_trainX = []
    n_trainY =[]
    zipxy = list(zip(trainX,trainY))
    random.shuffle(zipxy)
    for item1,item2 in zipxy:
        n_trainX.append(item1)
        n_trainY.append(item2)
    return n_trainX,n_trainY

def statistic(datalabel):
    kind2label = Counter(datalabel)
    for key,value in kind2label.items():
        yield "Label:{} Pro:{}".format(key,value/len(datalabel))

def classfication(trainX,trainY,testX,testY,model="randomforest",cv = 5,mode="cross_validate"):
    if model=="randomforest":
        clf = RandomForestClassifier(n_estimators=600, max_depth=None,
            min_samples_split=2, random_state=0)
    elif model =="bagging":
        clf = BaggingClassifier(DecisionTreeClassifier(random_state=0,max_depth=8))
    elif model=="boost":
        clf = AdaBoostClassifier(n_estimators=200,random_state=0)
    elif model=="vote":
        clf1 = RandomForestClassifier(n_estimators=500, max_depth=None,
            min_samples_split=2, random_state=0)
        clf2 = AdaBoostClassifier(n_estimators=100,random_state=0)
        clf3 = BaggingClassifier(DecisionTreeClassifier(random_state=0,max_depth=8))
        clf = VotingClassifier(estimators=[('rf', clf1), ('bg', clf3),('boost',clf2)], voting="soft")
    if mode=="predict":
        score = {}
        clf.fit(trainX,trainY)
        result = clf.predict(testX)
        p = precision_score(testY,result,average='macro')
        r = recall_score(testY,result,average='macro')
        f = f1_score(testY,result,average="macro")
        score["test_precision_macro"]=[p]
        score["test_recall_macro"]=[r]
        score["test_f1_macro"]=[f]
        return score
    elif mode=="cross_validate":
        X = []
        Y = []
        scoring = ["f1_macro", "precision_macro", "recall_macro"]
        X.extend(trainX)
        Y.extend(trainY)
        X.extend(testX)
        Y.extend(testY)
        score = cross_validate(clf, np.array(X),np.array(Y) , cv=cv, scoring=scoring, return_train_score=False)
        return score

def runpredict(trainX,trainY,testX,testY):
    score = classfication(trainX,trainY,testX,testY,model="randomforest",cv=5,mode="predict")
    f = score["test_f1_macro"]
    p = score["test_precision_macro"]
    r = score["test_recall_macro"]
    return  p,r,f
def runcross_validate(trainX,trainY,testX,testY):
    score = classfication(trainX,trainY,testX,testY,model="randomforest",cv=5,mode="cross_validate")
    f = np.mean(score["test_f1_macro"])
    p = np.mean(score["test_precision_macro"])
    r = np.mean(score["test_recall_macro"])
    return  p,r,f

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--disease2index",help="disease2index",required=True)
    parser.add_argument("--gene2index",help="gene2index",required=True)
    parser.add_argument("--intputall",help="intputall",required=True)
    parser.add_argument("--inputtrain",help="inputtrain",required=True)
    parser.add_argument("--outputA",help="outputA",required=True)
    args = parser.parse_args()
    disease2index = args.disease2index
    gene2index = args.gene2index
    intputall = args.intputall
    inputtrain = args.inputtrain
    outputA = args.outputA
    cfs=[]
    cps=[]
    crs=[]
    pfs=[]
    pps=[]
    prs=[]
    ratew = open("训练集比例.txt",'w')
    for i in tqdm(range(10)):
        trainX,trainY,testX,testY = load_data(intputall,inputtrain,outputA,disease2index,gene2index)
        ratew.write("Label distrubtion in train:")
        for item in statistic(trainY):
            ratew.write(item)

        ratew.write("Label distrubtion in test:")
        for item1 in statistic(testY):
            ratew.write(item1)


        cp,cr,cf = runcross_validate(trainX,trainY,testX,testY)
        pp,pr,pf = runpredict(trainX,trainY,testX,testY)
        cfs.append(cf)
        cps.append(cp)
        crs.append(cr)
        pfs.append(pf)
        pps.append(pp)
        prs.append(pr)
    
    print("-------------------------------cross_validate-------------------------------")
    print("P-value:{:.4f}±{:.4f}\tR-value:{:.4f}±{:.4f}\tF1-score:{:.4f}±{:.4f}".format(np.mean(cps),np.std(cps),np.mean(crs),np.std(crs),np.mean(cfs),np.std(cfs)))
    print("-------------------------------predict-------------------------------")
    print("P-value:{:.4f}±{:.4f}\tR-value:{:.4f}±{:.4f}\tF1-score:{:.4f}±{:.4f}".format(np.mean(pps),np.std(pps),np.mean(prs),np.std(prs),np.mean(pfs),np.std(pfs)))
if __name__=="__main__":
    main()