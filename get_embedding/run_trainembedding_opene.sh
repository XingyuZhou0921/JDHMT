#! usr/bin/bash
input_edges="/mnt/disk2/kyzhou/毕业设计第二部分/DataForWA/side_ediges.txt"
output_edges="/mnt/disk1/RESCALSIDE/GraphEmbedding/OpenNE/output/vector"
inputall="/mnt/disk2/kyzhou/毕业设计第二部分/triples/evidence.txt"
inputtrain="/mnt/disk2/kyzhou/毕业设计第二部分/triples/evidence_train.txt"
# inputtotal="/mnt/disk2/kyzhou/毕业设计第二部分/triples/gene_function_disease_total.txt"
# inputtrain="/mnt/disk2/kyzhou/毕业设计第二部分/triples/gene_function_disease_train.txt"
indexmapfile="/mnt/disk2/kyzhou/毕业设计第二部分/DataForWA/id2nodes.txt"
nodelabels="/mnt/disk2/kyzhou/毕业设计第二部分/DataForWA/nodelabels.txt"
representation_size=100
window_size=4
walk_length=30
lr=0.025
kstep=3 #for grarep
epochs=5 # for gf
weight_decay=0.0005
graph_format=edgelist
number_walks=20
data_mode="Unfilter"
arg=$1 #model name
echo
# echo "----------------------------------------------------------------------------------------------"
if [ $arg=="deepwalk" ];then
    echo Usage:
    echo python -m openne --method deepWalk --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --representation-size=$representation_size --window-size=$window_size --walk-length=$walk_length  --weighted --graph-format=$graph_format --lr=$lr
    echo "----------------------------------------------------------------------------------------------"
    python -m openne --method deepWalk --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --representation-size=$representation_size --window-size=$window_size --walk-length=$walk_length  --weighted --graph-format=$graph_format --lr=$lr

elif [ $arg=="line" ];then
    echo Usage:
    python -m openne --method line --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --window-size=$window_size --representation-size=$representation_size --order=3 --negative-ratio=5 --walk-length=$walk_length --weighted --graph-format=$graph_format --lr=$lr 
    echo "----------------------------------------------------------------------------------------------"
    python -m openne --method line --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --window-size=$window_size --representation-size=$representation_size --order=3 --negative-ratio=5 --walk-length=$walk_length --weighted --graph-format=$graph_format --lr=$lr 

elif [ $arg=="grarep" ];then
    echo Usage:
    python -m openne --method grarep --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --kstep=$kstep --representation-size=$representation_size --weighted --graph-format=$graph_format
    echo "----------------------------------------------------------------------------------------------"
    python -m openne --method grarep --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --kstep=$kstep --representation-size=$representation_size --weighted --graph-format=$graph_format 

elif [ $arg=="graphfactorization" ];then
    echo Usage:
    python -m openne --method gf --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --epochs=$epochs --weight-decay=$weight_decay --representation-size=$representation_size --weighted --graph-format=$graph_format --lr=$lr
    echo "----------------------------------------------------------------------------------------------"
    python -m openne --method gf --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --epochs=$epochs --weight-decay=$weight_decay --representation-size=$representation_size --weighted --graph-format=$graph_format --lr=$lr

elif [ $arg=="node2vec" ];then
    echo Usage:
    python -m openne --method node2vec --label-file=$nodelabels --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --graph-format=$graph_format --q=1 --p=0.25
    echo "----------------------------------------------------------------------------------------------"
    python -m openne --method node2vec --label-file=$nodelabels --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --graph-format=$graph_format --q=1 --p=0.25
    
elif [ $arg=="gcn" ];then
    echo Usage:
    python -m openne --method gcn --label-file=$nodelabels --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --graph-format=$graph_format --feature-file=data/cora/cora.features  --epochs=200 --clf-ratio=0.1
    echo "----------------------------------------------------------------------------------------------"
    # GCN method is a little different, 
    python -m openne --method gcn --label-file=$nodelabels --input=$input_edges --output=$output_edges"_"$arg"_"$data_mode --graph-format=$graph_format --feature-file=data/cora/cora.features  --epochs=200 --clf-ratio=0.1
fi

echo
echo "---------------Begin to evaluate the mode---------------" $arg

python classfication_OpenNE.py --inputtotal=$inputall --inputtrain=$inputtrain --vector=$output_edges"_"$arg"_"$data_mode --indexmapfile=$indexmapfile