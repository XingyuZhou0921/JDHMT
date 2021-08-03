#！ -*-encoding:utf-8-*-
from mapping import *

def get_gene_gene_relation(infile,outfile):
    wf = open(outfile,'w')
    rf = open(infile)
    wf.write("GeneA\tGeneB\n")
    rf.readline()
    for line in rf:
        contents = line.strip().split('\t')
        geneA = contents[1]
        geneB = contents[2]
        wf.write("{}\t{}\n".format(geneA,geneB))

def main():
    infile = "Gene_Gene/BIOGRID-ORGANISM-Homo_sapiens-3.5.179.tab2.txt"
    outfile = "Gene_Gene/gene_gene.txt"
    get_gene_gene_relation(infile,outfile)

if __name__=="__main__":
    main()
