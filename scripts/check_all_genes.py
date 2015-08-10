#!/usr/bin/env python

# This script takes 4 arguments as inputs: the data matrix, the cell labels (as determined by a clustering algorithm), names of all genes, names of potentially important genes. The script then attempts to figure out how well the clusters align with genes known to separate D1 cells from D2 cells. i.e. The script attempts to answer the question: Do the clusters correspond to D1 or D2 cells?

# Ideal results to confirm hypothesis that clusters correspond to D1/D2: Isl1, Drd1, Sfxn1, Nrxn1 are assigned to one cluster. Drd2, Penk, Sp9, Gpr52, Gpr88 are assigned to the other cluster. All p values are very small. 

import sys
import numpy as np
import scipy as sp
import scipy.stats
from tabulate import tabulate

def correct_nan(x):
    # if the value is NaN, change it to 1 (greater than any possibly log p value)
    if np.sum(np.isnan(x)) > 0:
        return 1
    else:
        return x

def chi_squared_test(count_0,count_1,labels,count_labels):
    # this function tests if the count distribution of samples between two clusters significantly deviates from no info
    p_obs = np.array([count_0,count_1])
    p_exp = np.array([count_labels[0],count_labels[1]]).astype(float)*(count_0+count_1)/np.sum(count_labels)
    test = sp.stats.chisquare(p_obs,f_exp=p_exp)
    return test[1]

def test_gene(clust_0,clust_1,thresh,labels,count_labels,gene):
    
    # number of samples in each cluster with expression less/greater than threshold
    count_0_less_thresh = np.sum(np.less(clust_0,thresh))
    count_0_greater_thresh = np.sum(np.greater(clust_0,thresh))
    count_1_less_thresh = np.sum(np.less(clust_1,thresh))
    count_1_greater_thresh = np.sum(np.greater(clust_1,thresh))

    # random assignly ties
    while count_0_less_thresh + count_0_greater_thresh < int(clust_0.shape[0]):
        if np.random.uniform(0,1,1) < count_labels[0]/np.sum(count_labels).astype(float):
            count_0_less_thresh += 1
        else:
            count_0_greater_thresh += 1
    while count_1_less_thresh +count_1_greater_thresh < int(clust_1.shape[0]):
        if np.random.uniform(0,1,1) < count_labels[0]/np.sum(count_labels).astype(float):
            count_1_less_thresh += 1
        else:
            count_1_greater_thresh += 1

    # run chi-squared test 
    p_using_greater = np.log(chi_squared_test(count_0_greater_thresh,count_1_greater_thresh,labels,count_labels))
    p_using_less = np.log(chi_squared_test(count_0_less_thresh,count_1_less_thresh,labels,count_labels))

    # correct for nans, just in case
    p_using_greater = correct_nan(p_using_greater)
    p_using_less = correct_nan(p_using_less)
    
    # use the smaller p-value
    test_results = []
    test_results.append(gene)
    if np.median(clust_0) > np.median(clust_1):
        test_results.append('0')
    else:
        test_results.append('1')
    test_results.append(np.min([p_using_greater,p_using_less]))

    return test_results     

def determine_clustering_quality(X,labels,uniq_labels,label_counts):
    # for each gene, extract the column of X that corresponds to expression of that gene
    save_data = []
    for i in range(0,X.shape[1]):

        gene_ind = i
        gene_exp = np.sort(X[:,gene_ind])
    
        # extract vector of gene expression for samples in each cluster
        clust_0 = X[labels==uniq_labels[0],gene_ind]
        clust_1 = X[labels!=uniq_labels[0],gene_ind]

        # test two thresholds: one assuming that the smaller cluster has less expression of the gene, and another assuming that the smaller cluster has more expression of the gene
        thresh_0 = np.mean([gene_exp[count_labels[0]-1],gene_exp[count_labels[0]]])
        thresh_1 = np.mean([gene_exp[count_labels[1]-1],gene_exp[count_labels[1]]])
        p_0 = test_gene(clust_0,clust_1,thresh_0,labels,count_labels,all_genes[i])
        p_1 = test_gene(clust_0,clust_1,thresh_1,labels,count_labels,all_genes[i])

        if p_0[2] > p_1[2]:
            add_data = p_1
        else:
            add_data = p_0

        save_data.append(add_data)
    return save_data
    # print '\n' + tabulate(save_data,headers=['Gene','Cluster Assignment','log(p-value)']) + '\n'

def extract_gene_stats(data,genes):
    qual = []
    for i in range(0,genes.shape[0]):
        gene_ind = int(np.where(all_genes==genes[i])[0])
        qual.append(data[gene_ind])
    return qual

# read in data matrix
X = np.loadtxt(sys.argv[1])

# read in labels
labels = np.loadtxt(sys.argv[2],dtype=str)
uniq_labels = np.unique(labels)
count_labels = np.array([np.sum(labels==uniq_labels[0]),np.sum(labels!=uniq_labels[0])])

# read in all genes
all_genes = np.loadtxt(sys.argv[3],dtype=str)

# read in D1/D2 genes
d1d2_genes = np.loadtxt('../data/genes_d1d2.txt',dtype=str)

# determine clustering quality for all genes!
qual = determine_clustering_quality(X,labels,uniq_labels,count_labels)

# extract clustering quality for d1d2 genes
qual_d1d2 = extract_gene_stats(qual,d1d2_genes)

if len(sys.argv) > 4:
    # read in important genes
    important_genes = np.loadtxt(sys.argv[4],dtype=str)
    qual_important_genes = extract_gene_stats(qual,important_genes)

# test random genes                                                                                                                                                             
rand_flag = sys.argv[5]
if rand_flag == 'htseq':
    rand_genes = np.loadtxt('/data/jessez/Gene_count_datasets/Genes/genes_rand_htseq.txt',dtype=str)
else:
    rand_genes = np.loadtxt('/data/jessez/Gene_count_datasets/Genes/genes_rand.txt',dtype=str)
qual_rand_genes = extract_gene_stats(qual,rand_genes)

# print stuff: also include L0 norm of cluster assignment to 000011111 or 111100000, mean of pvalues of D1/D2 genes
for i in range(0,len(qual_d1d2)-1):
    sys.stdout.write(qual_d1d2[i][0]+":")
sys.stdout.write(qual_d1d2[len(qual_d1d2)-1][0]+"\t")
for i in range(0,len(qual_d1d2)-1):
    sys.stdout.write(qual_d1d2[i][1])
sys.stdout.write(qual_d1d2[len(qual_d1d2)-1][1]+"\t")
for i in range(0,len(qual_d1d2)-1):
    sys.stdout.write(str(qual_d1d2[i][2])[0:7]+":")
sys.stdout.write(str(qual_d1d2[len(qual_d1d2)-1][2])[0:7]+"\t")

diff_000011111 = 0
diff_111100000 = 0
for i in range(0,4):
    if qual_d1d2[i][1] == '0':
        diff_111100000 += 1
    else:
        diff_000011111 += 1
for i in range(4,9):
    if qual_d1d2[i][1] == '0':
        diff_000011111 += 1
    else:
        diff_111100000 += 1
sys.stdout.write(str(np.min([diff_000011111,diff_111100000]))+"\t")

sum = 0
i = 1
for entry in qual_d1d2:
    sum += entry[2]
    i += 1
sys.stdout.write(str(sum/i)+"\t")

sum = 0
i = 1
for entry in qual_important_genes:
    sum += entry[2]
    i += 1
sys.stdout.write(str(sum/i)+"\t")

sum = 0
i = 1
for entry in qual:
    sum += entry[2]
    i += 1
sys.stdout.write(str(sum/i)+"\t")

sum = 0
i = 1
for entry in qual_rand_genes:
    sum += entry[2]
    i += 1
sys.stdout.write(str(sum/i))
