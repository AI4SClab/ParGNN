#!/bin/bash
#SBATCH -J products
#SBATCH -p normal # debug # normal
#SBATCH -N 4
#SBATCH --gres=dcu:4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --exclusive
#SBATCH --mem=90G
#SBATCH --time=10:00:00
#SBATCH --output=./logs/train_exp_%j.out

mkdir -p ./logs
# source environment for the cluster
source ./xd_start.sh
which python

export MASTER_PORT=25875
export MASTER_ADDR=$(scontrol show hostname ${SLURM_NODELIST} | head -n 1)
export NCCL_SHM_DISABLE=1
export NCCL_BLOCKING_WAIT=1
export NCCL_ASYNC_ERROR_HANDLING=1
export NCCL_P2P_LEVEL=PXB # PHB
export GLOO_DEBUG=INFO
# echo $(ib0)
# echo $(eth0)
export GLOO_SOCKET_IFNAME=ib0
export NCCL_SOCKET_IFNAME=ib0

echo start on $(date)
echo "SLURM_JOB_ID: $SLURM_JOB_ID" 


# Training set up
cd ../src/

num_gpu=${SLURM_NTASKS}
dulpicator=4
backend="gloo"

epochs=50
total_parts=$((num_gpu * dulpicator))
seed=2026
data_path="../data/DistData"


# Define datasets to run: "dataset_name:num_classes"
datasets=(
    "ogbn-products:47"
    "ogbn-proteins:112"
    "ogbn-arxiv:40"
    "yelp:100"
    "reddit:41"
)

for item in "${datasets[@]}"; do
    # Split dataset name and num_class
    dataset="${item%%:*}"
    num_class="${item##*:}"
    
    echo "=========================================="
    echo "Training on dataset: $dataset (num_class=$num_class)"
    echo "=========================================="
    gcn_app_0_base="python main_train.py \
        --epochs $epochs \
        --hidden 256 \
        --model gcn \
        --num_layers 3 \
        --total_parts $total_parts \
        --seed $seed \
        --data_path $data_path \
        --backend $backend \
        --dataset $dataset \
        --num_class $num_class \
        --use_pipeline 0 \
        --reparter base"



    gcn_app_1_base="python main_train.py \
        --epochs $epochs \
        --hidden 256 \
        --model gcn \
        --num_layers 3 \
        --total_parts $total_parts \
        --seed $seed \
        --data_path $data_path \
        --backend $backend \
        --dataset $dataset \
        --num_class $num_class \
        --use_pipeline 1 \
        --reparter base"
    
    gcn_app_1_adapt="python main_train.py \
        --epochs $epochs \
        --hidden 256 \
        --model gcn \
        --num_layers 3 \
        --total_parts $total_parts \
        --seed $seed \
        --data_path $data_path \
        --backend $backend \
        --dataset $dataset \
        --num_class $num_class \
        --use_pipeline 1 \
        --reparter adapt"



    gat_app_0_base="python main_train.py \
        --epochs $epochs \
        --hidden 256 \
        --model gat \
        --num_layers 3 \
        --total_parts $total_parts \
        --seed $seed \
        --data_path $data_path \
        --backend $backend \
        --dataset $dataset \
        --num_class $num_class \
        --use_pipeline 0 \
        --reparter base"

    
    gat_app_1_base="python main_train.py \
        --epochs $epochs \
        --hidden 256 \
        --model gat \
        --num_layers 3 \
        --total_parts $total_parts \
        --seed $seed \
        --data_path $data_path \
        --backend $backend \
        --dataset $dataset \
        --num_class $num_class \
        --use_pipeline 1 \
        --reparter base"


    
    gat_app_1_adapt="python main_train.py \
        --epochs $epochs \
        --hidden 256 \
        --model gat \
        --num_layers 3 \
        --total_parts $total_parts \
        --seed $seed \
        --data_path $data_path \
        --backend $backend \
        --dataset $dataset \
        --num_class $num_class \
        --use_pipeline 1 \
        --reparter adapt"

    

    
    mpirun -n $num_gpu $gcn_app_0_base
    echo "Time: $(date)"
    mpirun -n $num_gpu $gcn_app_1_adapt
    echo "Time: $(date)" 

    mpirun -n $num_gpu $gat_app_0_base
    echo "Time: $(date)"
    mpirun -n $num_gpu $gat_app_1_adapt
    echo "Time: $(date)"
    echo "Finished training on $dataset"
    echo "Time: $(date)"
    echo ""
done