#! /usr/bin/bash
A=0.01
R=0.01
V=0.01
S=0.001
init_method="random"
fit_method=False #whether we need to calculate the 
checkpoint="1000,1100" #Calculate fit
maxiter=10 # for other models
inputall="/mnt/disk2/kyzhou/毕业设计第二部分/triples/evidence.txt"
inputtrain="/mnt/disk2/kyzhou/毕业设计第二部分/triples/evidence_train.txt"
# inputall="/mnt/disk2/kyzhou/毕业设计第二部分/triples/gene_function_disease_total.txt"
# inputtrain="/mnt/disk2/kyzhou/毕业设计第二部分/triples/gene_function_disease_train.txt"

version=5 #use rescalv3 or rescalv4

echo $1
echo "2019/12/31"

outputA="output/A9.npy"
outputR="output/R9.npy"
outputX="output/X9.npy"
disease2indexpath="output/disease2index9.npy"
gene2indexpath="output/gene2index9.npy"

echo "----------------------------------Paramater Setting----------------------------------"
echo "Using lambdaA" $A "using lambdaR" $R "using lambdaV" $V 'using lambdaS' $S
echo 'using Maxiter' $maxiter "init_method" $init_method "fit_method" $fit_method "checkpoint" $checkpoint 
python runonlypca.py --lambdaA=$A --lambdaR=$R --lambdaV=$V --lambdaS=$S --init_method=$init_method --fit_method=$fit_method --checkpoint=$checkpoint --maxiter=$maxiter --outputA=$outputA --outputR=$outputR --outputX=$outputX --disease2indexpath=$disease2indexpath --gene2indexpath=$gene2indexpath --inputall=$inputall --inputtrain=$inputtrain --model="PCA"

echo
echo "--------------------------Begin to evaluate the mode-------------------------" $arg
python classfication_Rescal.py --disease2index=$disease2indexpath --gene2index=$gene2indexpath --intputall=$inputall --inputtrain=$inputtrain --outputA=$outputA
echo
