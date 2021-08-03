#! -*- encoding:utf-8 -*-
from mapping import *

def get_disease_disease(inputpath,outputpath,diseaseontolgy):
    do2mesh = disease_DO2mesh(diseaseontolgy)
    count=0
    with open(inputpath) as rf:
        with open(outputpath,'w') as wf:
            rf.readline()
            wf.write("DiseaseA\tDiseaseB\n")
            for line in rf:
                contents = line.strip().split('\t')
                disease1 = do2mesh.get(contents[0],-10)
                disease = do2mesh.get(contents[1],-10)
                if disease!=-10 and disease1!=-10:
                    wf.write("{}\t{}\n".format(disease1,disease))
                else:
                    count+=1
    print("There are {} we can't map!".format(count))
    
def main():
    diseaseontolgy="vocabs/diseaseontolgy.obo"
    inputpath = "Disease_Disease/DD-Miner_miner-disease-disease.tsv"
    outputpath = "Disease_Disease/disease_disease.txt"
    get_disease_disease(inputpath,outputpath,diseaseontolgy)

if __name__=="__main__":
    main()