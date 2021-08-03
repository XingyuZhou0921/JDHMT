#! -*- encoding:utf-8 -*-
import re
import codecs
def get_gene_disease(fileinput,fileoutput):
    wf = codecs.open(fileoutput,'w',encoding='utf-8')
    wf.write("GeneID\tDiseaseID\n")
    with codecs.open(fileinput,encoding='utf-8') as rf:
        for line in rf:
            if line.startswith("#"):
                pass
            else:
                line = re.sub('"(.*?)"',"REPLACE",line)
                contnts = line.strip().split(',')
                geneid = contnts[1]
                diseaseid = contnts[3]
                wf.write("{}\t{}\n".format(geneid,diseaseid))


def main():
    fileinput="gene-disease/CTD_genes_diseases.csv"
    fileroutput = "gene-disease/gene_disease.txt"

    get_gene_disease(fileinput,fileroutput)

if __name__=="__main__":
    main()
