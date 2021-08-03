#! /usr/bin/bash
A=0.0001
R=0.0001
V=0.0001
S=0.001
init_method="random"
fit_method=False #whether we need to calculate the 
checkpoint="1000,1100" #Calculate fit
maxiter=300
inputall="/mnt/disk2/kyzhou/毕业设计第二部分/triples/evidence.txt"
inputtrain="/mnt/disk2/kyzhou/毕业设计第二部分/triples/evidence_train.txt"
# inputall="/mnt/disk2/kyzhou/毕业设计第二部分/triples/gene_function_disease_total.txt"
# inputtrain="/mnt/disk2/kyzhou/毕业设计第二部分/triples/gene_function_disease_train.txt"

version=4 #use rescalv3 or rescalv4

echo "2019/1/4"

outputA="output/A11.npy"
outputR="output/R11.npy"
outputX="output/X11.npy"
disease2indexpath="output/disease2index11.npy"
gene2indexpath="output/gene2index11.npy"

# echo "----------------------------------Paramater Setting----------------------------------"
echo "Using lambdaA" $A "using lambdaR" $R "using lambdaV" $V 'using lambdaS' $S
echo 'using Maxiter' $maxiter "init_method" $init_method "fit_method" $fit_method "checkpoint" $checkpoint 
python runrescalv$version.py --lambdaA=$A --lambdaR=$R --lambdaV=$V --lambdaS=$S --init_method=$init_method --fit_method=$fit_method --checkpoint=$checkpoint --maxiter=$maxiter --outputA=$outputA --outputR=$outputR --outputX=$outputX --disease2indexpath=$disease2indexpath --gene2indexpath=$gene2indexpath --inputall=$inputall --inputtrain=$inputtrain

echo
echo "--------------------------Begin to evaluate the mode-------------------------" $arg
python classfication_Rescal.py --disease2index=$disease2indexpath --gene2index=$gene2indexpath --intputall=$inputall --inputtrain=$inputtrain --outputA=$outputA
echo "Paramater:Maxiterator" $maxiter
echo

