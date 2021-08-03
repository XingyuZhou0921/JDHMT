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

echo $1
echo "2019/12/30"
if [[ $1 = 'maxiter' ]];then
    for maxiter in {50..550..50}
        do
            outputA="output/A4.npy"
            outputR="output/R4.npy"
            outputX="output/X4.npy"
            disease2indexpath="output/disease2index4.npy"
            gene2indexpath="output/gene2index4.npy"

            echo "----------------------------------Paramater Setting----------------------------------"
            echo "Using lambdaA" $A "using lambdaR" $R "using lambdaV" $V 'using lambdaS' $S
            echo 'using Maxiter' $maxiter "init_method" $init_method "fit_method" $fit_method "checkpoint" $checkpoint 
            python runrescalv$version.py --lambdaA=$A --lambdaR=$R --lambdaV=$V --lambdaS=$S --init_method=$init_method --fit_method=$fit_method --checkpoint=$checkpoint --maxiter=$maxiter --outputA=$outputA --outputR=$outputR --outputX=$outputX --disease2indexpath=$disease2indexpath --gene2indexpath=$gene2indexpath --inputall=$inputall --inputtrain=$inputtrain

            echo
            echo "--------------------------Begin to evaluate the mode-------------------------" $arg
            python classfication_Rescal.py --disease2index=$disease2indexpath --gene2index=$gene2indexpath --intputall=$inputall --inputtrain=$inputtrain --outputA=$outputA
            echo "Paramater:Maxiterator" $maxiter
            echo
        done
elif [[ $1 = 'ARV' ]];then
    for V in 0.00001 0.0001 0.001 0.01 0.1 1
        do  
            maxiter=200
            outputA="output/A5.npy"
            outputR="output/R5.npy"
            outputX="output/X5.npy"
            disease2indexpath="output/disease2index5.npy"
            gene2indexpath="output/gene2index5.npy"

            echo "----------------------------------Paramater Setting----------------------------------"
            echo "Using lambdaA" $V "using lambdaR" $V "using lambdaV" $V 'using lambdaS' $S
            echo 'using Maxiter' $maxiter "init_method" $init_method "fit_method" $fit_method "checkpoint" $checkpoint 
            python runrescalv$version.py --lambdaA=$V --lambdaR=$V --lambdaV=$V --lambdaS=$S --init_method=$init_method --fit_method=$fit_method --checkpoint=$checkpoint --maxiter=$maxiter --outputA=$outputA --outputR=$outputR --outputX=$outputX --disease2indexpath=$disease2indexpath --gene2indexpath=$gene2indexpath --inputall=$inputall --inputtrain=$inputtrain

            echo
            echo "--------------------------Begin to evaluate the mode-------------------------" $arg
            python classfication_Rescal.py --disease2index=$disease2indexpath --gene2index=$gene2indexpath --intputall=$inputall --inputtrain=$inputtrain --outputA=$outputA
            echo "Paramater:lambdaA=lambdaR=lambdaV" $lambdaA
            echo
        done
elif [[ $1 = 'SS' ]];then # paramater lamadaS
    for S in 0.00001 0.0001 0.001 0.01 0.1 1
        do  
            maxiter=200
            outputA="output/A6.npy"
            outputR="output/R6.npy"
            outputX="output/X6.npy"
            disease2indexpath="output/disease2index6.npy"
            gene2indexpath="output/gene2index6.npy"

            echo "----------------------------------Paramater Setting----------------------------------"
            echo "Using lambdaA" $A "using lambdaR" $R "using lambdaV" $V 'using lambdaS' $S
            echo 'using Maxiter' $maxiter "init_method" $init_method "fit_method" $fit_method "checkpoint" $checkpoint 
            python runrescalv$version.py --lambdaA=$A --lambdaR=$R --lambdaV=$V --lambdaS=$S --init_method=$init_method --fit_method=$fit_method --checkpoint=$checkpoint --maxiter=$maxiter --outputA=$outputA --outputR=$outputR --outputX=$outputX --disease2indexpath=$disease2indexpath --gene2indexpath=$gene2indexpath --inputall=$inputall --inputtrain=$inputtrain

            echo
            echo "--------------------------Begin to evaluate the mode-------------------------" $arg
            python classfication_Rescal.py --disease2index=$disease2indexpath --gene2index=$gene2indexpath --intputall=$inputall --inputtrain=$inputtrain --outputA=$outputA
            echo "Paramater:lambdaS" $lambdaS
            echo
        done
else 
    echo "Input Error!"
fi
