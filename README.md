# JDHMT
JDHMT: High-quality Gene/Disease Embedding in A Multi-relational Heterogeneous Graph After A Joint Matrix/tensor Decomposition

## Tested environment
**Python>=3.5**

## Basic Usage
__1.Get Primary data(file:json)__  
Step(1): run ./get_primarydata/1 to 10 .py(which likes __1_chemicalrelation.py__)  
input = './CTD_chemicals_disease.csv'  
output = './chemical_disease.txt'  
Step(2): run ./get_primarydata/build_graph_format.py  
input = './gene_function_disease_train.txt' ,'side_ediges.txt','nodelabels.txt' 
output = './DataForWA/id2nodes.txt'  
Step(2): run ./get_primarydata/build_w_json.py  
input = './(diseasename/genename).txt' 
output = './chemical_gene.txt' 

