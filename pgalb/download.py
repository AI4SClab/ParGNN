import dgl
import numpy as np
import torch, json
import os
import os.path as osp
import torch.multiprocessing as mp
from typing import List, Dict
import time
from ogb.nodeproppred import DglNodePropPredDataset
def call_back(res):
    return None

def err_call_back(err):
    print(f'!!! ~ error: {str(err)}')

def load_check_download(
    name:str = "ogbn-arxiv",
    dir_path:str = "/work/gujy/data",
    save_path:str = "/work/gujy/data/DistData"
):
    """
    Load dataset for check and download datasets
    """
    if name in ['ogbn-arxiv','ogbn-products']:
        dataset = DglNodePropPredDataset(name=name,root=dir_path)
        graph : dgl.DGLGraph = dataset[0][0]

        feat = graph.ndata['feat']
        label : torch.Tensor = dataset[0][1].reshape(graph.num_nodes(),-1) # type: ignore
        split_idx  = dataset.get_idx_split()
        mask = torch.ones(size=(graph.num_nodes(),), dtype=torch.int32).mul_(-1)
        mask[split_idx["train"]] = 0 # type: ignore
        mask[split_idx["valid"]] = 1 # type: ignore
        mask[split_idx["test"]] = 2  # type: ignore
    elif name in ['ogbn-proteins']:
        dataset = DglNodePropPredDataset(name='ogbn-proteins',root=dir_path)
        # print(dataset[0][0])
        graph : dgl.DGLGraph = dataset[0][0]
        label : torch.Tensor = dataset[0][1].reshape(graph.num_nodes(),-1) # type: ignore
        A = graph.adj()
        # A = dgl.sparse.val_like(A,graph.edata['feat'])
        # feat = A.smean(1)
        split_idx  = dataset.get_idx_split()
        mask = torch.ones(size=(graph.num_nodes(),), dtype=torch.int32).mul_(-1)
        mask[split_idx["train"]] = 0 # type: ignore
        mask[split_idx["valid"]] = 1 # type: ignore
        mask[split_idx["test"]] = 2  # type: ignore
    elif name in ['yelp','reddit', 'flickr']:
        from dgl.data import YelpDataset, RedditDataset
        dataset = YelpDataset(raw_dir=dir_path) if name == 'yelp' else RedditDataset(raw_dir=dir_path)
        if name == 'flickr':
            from dgl.data import FlickrDataset
            dataset = FlickrDataset(raw_dir=dir_path)
        
        dataset.num_classes
        graph = dataset[0]
        # print(graph)
        # get node feature
        feat = graph.ndata['feat']
        # get node labels
        label = graph.ndata['label']
        # get data split
        train_mask = graph.ndata['train_mask'].to(torch.bool)
        val_mask = graph.ndata['val_mask'].to(torch.bool)
        test_mask = graph.ndata['test_mask'].to(torch.bool)

        mask = torch.ones(size=(graph.num_nodes(),), dtype=torch.int32).mul_(-1)
        mask[train_mask] = 0
        mask[val_mask] = 1
        mask[test_mask] = 2

    elif 'igb' in name:
        _ , dataset_size = name.split('-')
        from igb.dataloader import IGB260M
        dir_path = osp.join(dir_path,"IGB260M")
        ## need to download first
        dataset = IGB260M(root=dir_path, size=dataset_size, in_memory=1, classes=19, synthetic=0)
        feat = torch.from_numpy(dataset.paper_feat)
        node_edges = torch.from_numpy(dataset.paper_edge)
        label = torch.from_numpy(dataset.paper_label).to(torch.long)

        graph = dgl.graph((node_edges[:, 0],node_edges[:, 1]), num_nodes=feat.shape[0])

        mask = torch.ones(size=(graph.num_nodes(),), dtype=torch.int32).mul_(-1)
        n_nodes = feat.shape[0]
        n_train = int(n_nodes * 0.6)
        n_val   = int(n_nodes * 0.2)

        mask[:n_train] = 0
        mask[n_train:n_train + n_val] = 1
        mask[n_train + n_val:] = 2 

    else:
        
        print(f"Data: {name} is not support now, please check it!")
        raise 
    return 

if __name__ == '__main__':

    datasets : List[str] = ['yelp','reddit',"ogbn-arxiv",'ogbn-products','ogbn-proteins']
    dir_path : str = "../data"
    save_path: str = "../data/DistData"
    for name in datasets:
        load_check_download(name=name, dir_path= dir_path, save_path = save_path )