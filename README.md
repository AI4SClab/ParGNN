# ParGNN: A Scalable Graph Neural Network Training Framework on multi-GPUs

ParGNN is accepted by DAC 2025.

ParGNN, an efficient full-batch training system for GNNs, which adopts a profiler-guided adaptive load balancing partition method(PGALB) and a subgraph pipeline algorithm to overlap communication and computation.

### 1. Clone this project and Setup environment

### 2. Install PGALB and partition graph data

###### 2.1. setup pgalb

```shell
# [update] Fixed the miss of source code for GKlib and METIS
cd ParGNN # make sure you are in the ParGNN directory
export WORK_DIR=$(pwd)

cd $WORK_DIR/pgalb/csrc/GKlib
rm -rf build/
make config prefix=$WORK_DIR/pgalb/csrc/local
make install

cd $WORK_DIR/pgalb/csrc/METIS
rm -rf build/
make config prefix=$WORK_DIR/pgalb/csrc/local gklib_path=$WORK_DIR/pgalb/csrc/local i64=1
make install

cd $WORK_DIR/pgalb
python setup.py build_ext --inplace  ## install pgalb
python python test_adapt.py   ## check the install of C extension

python download.py
```

###### 2.2. parition and reparition

* ParGNN use two stage partition to deal with the load unbalance question, and provide the 3 functions:
    * graph_partition_dgl_metis:  load graph from local path(or download online automately), and create DGL-format graph. Then run the intial patition by metis.
    * graph_eval: profile the subgraphs from the inital parition.
    * mapping: map the subgraphs to a cluster for running in the same GPUs.
* we provide two scripts to test the partition and repartition process, just run:
    ```shell
    python graph_partition.py
    python repart.py
    ```
* data used in the evaluation on paper: (script will download them automately when they are needed)
    * ogb graph dataset: ogbn-products, ogbn-proteins from https://ogb.stanford.edu/docs/nodeprop/
    * yelp dataset: https://www.dgl.ai/dgl_docs/generated/dgl.data.YelpDataset.html#dgl.data.YelpDataset
    * reddit dataset: https://www.dgl.ai/dgl_docs/generated/dgl.data.RedditDataset.html

### 3. Run distributed GNN trainning

* The scripts dirctionary has the example scripts to run ParGNN. 
    ```shell
    cd scripts
    sh train_all.sh
    ```
### 4. [update] Support for slurm scheduler for GPU cluster on Multi-Node training

* The slurm scripts are in the slurm_scripts directory.
    * `slurm_scripts/partition_and_repart.slurm`: partition and repartition script for slurm scheduler
    * `slurm_scripts/train.slurm`: train script for slurm scheduler

### 5. Citations

To cite this project, you can use the following BibTex citation.

[DAC 2025](https://ieeexplore.ieee.org/document/11133102 "dac2025")
```
@INPROCEEDINGS{11133102,
  author={Gu, Junyu and Li, Shunde and Cao, Rongqiang and Wang, Jue and Wang, Zijian and Liang, Zhiqiang and Liu, Fang and Li, Shigang and Zhou, Chunbao and Wang, Yangang and Chi, Xuebin},
  booktitle={2025 62nd ACM/IEEE Design Automation Conference (DAC)}, 
  title={ParGNN: A Scalable Graph Neural Network Training Framework on multi-GPUs}, 
  year={2025},
  volume={},
  number={},
  pages={1-7},
  keywords={Training;Accuracy;Design automation;Pipelines;Graphics processing units;Load management;Graph neural networks;Partitioning algorithms;Faces;Convergence;Graph neural network;Full-batch distributed training;Load balancing;Computation and communication overlapping},
  doi={10.1109/DAC63849.2025.11133102}}
```

[PPoPP 2024 Poster](https://dl.acm.org/doi/abs/10.1145/3627535.3638488 "ppopp2024")
```
@inproceedings{10.1145/3627535.3638488,
author = {Li, Shunde and Gu, Junyu and Wang, Jue and Yao, Tiechui and Liang, Zhiqiang and Shi, Yumeng and Li, Shigang and Xi, Weiting and Li, Shushen and Zhou, Chunbao and Wang, Yangang and Chi, Xuebin},
title = {POSTER: ParGNN: Efficient Training for Large-Scale Graph Neural Network on GPU Clusters},
year = {2024},
isbn = {9798400704352},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3627535.3638488},
doi = {10.1145/3627535.3638488},
pages = {469–471},
numpages = {3},
keywords = {graph neural network, load balancing, data transfer hiding, distributed training},
location = {Edinburgh, United Kingdom},
series = {PPoPP '24}
}
```
































































*
