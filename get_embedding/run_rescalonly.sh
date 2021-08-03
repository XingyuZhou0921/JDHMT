#! /usr/bin/bash
A=0.01
R=0.01
V=0.01
S=0.001
init_method="random"
fit_method=False #whether we need to calculate the 
checkpoint="1000,1100" #Calculate fit

inputall="/mnt/disk2/kyzhou/毕业设计第二部分/triples/evidence.txt"
inputtrain="/mnt/disk2/kyzhou/毕业设计第二部分/triples/evidence_train.txt"
# inputall="/mnt/disk2/kyzhou/毕业设计第二部分/triples/gene_function_disease_total.txt"
# inputtrain="/mnt/disk2/kyzhou/毕业设计第二部分/triples/gene_function_disease_train.txt"

version=5 #use rescalv3 or rescalv4
maxiter=300

outputA="output/A27.npy"
outputR="output/R27.npy"
outputX="output/X27.npy"
disease2indexpath="output/disease2index27.npy"
gene2indexpath="output/gene2index27.npy"

# 注意这边要调整 count和0/1模式


echo "----------------------------------Paramater Setting----------------------------------"
echo "Using lambdaA" $A "using lambdaR" $R "using lambdaV" $V 'using lambdaS' $S
echo 'using Maxiter' $maxiter "init_method" $init_method "fit_method" $fit_method "checkpoint" $checkpoint 
python runonlyrescal.py --lambdaA=$A --lambdaR=$R --lambdaV=$V --lambdaS=$S --init_method=$init_method --fit_method=$fit_method --checkpoint=$checkpoint --maxiter=$maxiter --outputA=$outputA --outputR=$outputR --outputX=$outputX --disease2indexpath=$disease2indexpath --gene2indexpath=$gene2indexpath --inputall=$inputall --inputtrain=$inputtrain --model="RESCAL"

echo
echo "--------------------------Begin to evaluate the mode-------------------------" $arg
python classfication_Rescal.py --disease2index=$disease2indexpath --gene2index=$gene2indexpath --intputall=$inputall --inputtrain=$inputtrain --outputA=$outputA
echo "Paramater:Maxiterator" $maxiter
echo
       
