# gene and disease relations with different bioconcepts

## Chemical–disease associations 
Fields:
* 1    ChemicalID (MeSH identifier)
* 2    DiseaseID (MeSH identifier)
* 3    OmimIDs ('|'-delimited list)
* 4    PubMedIDs ('|'-delimited list)

## Chemical–gene interactions
Fields:
* 1    ChemicalID (MeSH identifier)
* 2    GeneID (NCBI Gene identifier)
* 3    InteractionActions ('|'-delimited list)
* 4    PubMedIDs ('|'-delimited list)

## Pathway-Gene associations
Fields:
* 1    PathwayID (KEGG or REACTOME identifier)
* 2    GeneID (NCBI Gene identifier)

## Pathway-Disease associations		
Fields:
* 1    PathwayID (KEGG or REACTOME identifier)
* 2    DiseaseID (MeSH or OMIM identifier)


## Phenotype-Disease associations		
Fields:
* 1	 PhenotypeID (HPO)
* 2	 DiseaseID (MeSH or OMIM identifier)	
	
	
## Phenotype-Gene associations		
Fields:
* 1	 PathwayID (HPO)
* 2	 GeneID (NCBI Gene identifier)		
	
	
## Mutation-Disease associations		
Fields:
* 1	 VariationID （Clinvar variation identifier）              
* 2	 PhenotypeIDs  (MedGen,OMIM,HPO,SNOMED CT,Orphanet)
(here we could only get the Mutation-phenotypeIDs. phenotypeIDs not equal disease)        


## Mutation-Gene associations		
Fields:
* 1	 VariationID (Clinvar variation identifier)
* 2	 GeneID (NCBI Gene identifier)	



# Triples for knowledge inference

## Triples:  "Gene-Function_change-Disease"

It is a fact that functional changes occur after gene mutations may lead to disease. 
By mining the existing triples (Gene-function_change-disease) and reasoning out new relations,
that will provide precise guidance for drug development.[1]


Current:
* we have designed a pipeline that could help obtain 'gene-function_change-disease' triples from PubMed abstracts.
* We want to design relational reasoning algorithms around tensor decomposition to acquire unknown hidden knowledge. 

Some thoughts：
* Using and Improving Tensor Decomposition Algorithms for Relational Reasoning(such as RESCAL).
* Combine Tensor Decomposition with neural network or traditional machine learning methods.

## Triples
Fields:
* 1   PMID
* 2   Gene                 (should map to id)
* 3   Function-change     （GOF/LOF/COM     GOF means loss of function, GOF means Gain of function, COM means this gene's different mutations can cause LOF or GOF）
* 4   Disease             （should map to id）


Only contain two differen disease: bladder cancer(BLCA)


[1] Li, Yongsheng, et al. "Gain-of-Function Mutations: An Emerging Advantage for Cancer Biology." Trends in Biochemical Sciences. 



















