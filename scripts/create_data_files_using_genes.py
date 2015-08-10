#!/bin/bash/python

import sys
import os
import numpy as np

def get_gene_indices(all_genes,genes):
    ind = np.zeros(len(all_genes))
    for i in range(0,len(genes)):
        ind += (all_genes == genes[i])
    return ind

rand_genes = np.loadtxt("/data/jessez/Gene_count_datasets/Genes/genes_rand.txt",dtype=str)
rand_genes_htseq = np.loadtxt("/data/jessez/Gene_count_datasets/Genes/genes_rand_htseq.txt",dtype=str)
all_genes = np.loadtxt("/data/jessez/Gene_count_datasets/Genes/genes.txt",dtype=str)
all_genes_htseq = np.loadtxt("/data/jessez/Gene_count_datasets/Genes/genes_htseq.txt",dtype=str)

rand_genes_ind = np.where(get_gene_indices(all_genes,rand_genes))
rand_genes_ind_htseq = np.where(get_gene_indices(all_genes_htseq,rand_genes_htseq))

data_dir = "/data/jessez/Gene_count_datasets/Datasets"
for i in os.listdir(data_dir):
    for j in os.listdir(data_dir+"/"+i):
        if "d1d2" not in j and "Meisam" not in j:
            X = np.loadtxt(data_dir+"/"+i+"/"+j)
            if "htseq" in j:
                X = np.squeeze(X[:,rand_genes_ind_htseq])
            else:
                X = np.squeeze(X[:,rand_genes_ind])
            np.savetxt(data_dir+"/"+i+"/"+j+"_RandGenes",X,delimiter='\t')
