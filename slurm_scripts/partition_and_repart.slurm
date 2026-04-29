#!/bin/bash
#SBATCH -J partition_and_repart
#SBATCH -p normal # debug # normal
#SBATCH -N 1
#SBATCH --gres=dcu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=90G
#SBATCH --output=./logs/part_repart_%j.out

# source environment for the cluster
source ./xd_start


cd ../pgalb

echo "Current directory: $(pwd)"
echo "Starting graph partitioning at date: $(date)"
echo "----------------------------------------"
python graph_partition.py
echo "----------------------------------------"
echo "Starting graph repartitioning at date: $(date)"
echo "----------------------------------------"
python repart.py
echo "----------------------------------------"
echo "Finished at date: $(date)"

