#!/usr/bin/env python

# This script takes 4 arguments as inputs: the data matrix, the cell labels (as determined by a clustering algorithm), names of all genes, names of potentially important genes. The script then attempts to figure out how well the clusters align with genes known to separate D1 cells from D2 cells. i.e. The script attempts to answer the question: Do the clusters correspond to D1 or D2 cells?

# Ideal results to confirm hypothesis that clusters correspond to D1/D2: Isl1, Drd1, Sfxn1, Nrxn1 are assigned to one cluster. Drd2, Penk, Sp9, Gpr52, Gpr88 are assigned to the other cluster. All p values are very small. 

import sys
import numpy as np
import scipy as sp
import scipy.stats
from tabulate import tabulate

def chi_squared_test(count_0,count_1,labels,count_labels):
    # this function tests if the count distribution of samples between two clusters significantly deviates from no info
    p_obs = np.array([count_0,count_1])
    p_exp = np.array([count_labels[0],count_labels[1]]).astype(float)*(count_0+count_1)/np.sum(count_labels)
    test = sp.stats.chisquare(p_obs,f_exp=p_exp)
    # print p_obs
    # print p_exp
    # print test[1]
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
    # print '  with greater'
    p_using_greater = np.log(chi_squared_test(count_0_greater_thresh,count_1_greater_thresh,labels,count_labels))
    # print '  with less'
    p_using_less = np.log(chi_squared_test(count_0_less_thresh,count_1_less_thresh,labels,count_labels))
    
    # use the smaller p-value
    test_results = []
    test_results.append(gene)
    if np.mean(clust_0) > np.mean(clust_1):
        test_results.append('0')
    else:
        test_results.append('1')
    test_results.append(np.min([p_using_greater,p_using_less]))

    return test_results     

def determine_clustering_quality(X,labels,uniq_labels,label_counts,genes):
    # for each gene, extract the column of X that corresponds to expression of that gene
    save_data = []
    for i in range(0,int(genes.shape[0])):

        # print '\n'

        gene_ind = int(np.where(all_genes==genes[i])[0])
        gene_exp = np.sort(X[:,gene_ind])
    
        # extract vector of gene expression for samples in each cluster
        clust_0 = X[labels==uniq_labels[0],gene_ind]
        clust_1 = X[labels!=uniq_labels[0],gene_ind]

        # test two thresholds: one assuming that the smaller cluster has less expression of the gene, and another assuming that the smaller cluster has more expression of the gene
        thresh_0 = np.mean([gene_exp[count_labels[0]-1],gene_exp[count_labels[0]]])
        thresh_1 = np.mean([gene_exp[count_labels[1]-1],gene_exp[count_labels[1]]])

        # print 'testing thresh 1'
        p_0 = test_gene(clust_0,clust_1,thresh_0,labels,count_labels,genes[i])
        # print 'testing thresh 2'
        p_1 = test_gene(clust_0,clust_1,thresh_1,labels,count_labels,genes[i])

        if p_0[2] > p_1[2]:
            add_data = p_1
        else:
            add_data = p_0

        save_data.append(add_data)
    return save_data
    # print '\n' + tabulate(save_data,headers=['Gene','Cluster Assignment','log(p-value)']) + '\n'


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

# determine clustering quality!
qual = determine_clustering_quality(X,labels,uniq_labels,count_labels,d1d2_genes)

if len(sys.argv) > 4:
    # read in important genes
    important_genes = np.loadtxt(sys.argv[4],dtype=str)
    qual_important_genes = determine_clustering_quality(X,labels,uniq_labels,count_labels,important_genes[np.random.permutation(99)[0:9]])

# print stuff version 1: gene names, cluster assignments, pvalues (and mean of pvalues of important genes)
for i in range(0,len(qual)-1):
    sys.stdout.write(qual[i][0]+":")
sys.stdout.write(qual[len(qual)-1][0]+"\t")
for i in range(0,len(qual)-1):
    sys.stdout.write(qual[i][1])
sys.stdout.write(qual[len(qual)-1][1]+"\t")
for i in range(0,len(qual)-1):
    sys.stdout.write(str(qual[i][2])[0:7]+":")
sys.stdout.write(str(qual[len(qual)-1][2])[0:7]+"\t")

# sum = 0
# i = 1
# for entry in qual_important_genes:
#     sum += entry[2]
#     i += 1
# sys.stdout.write(str(sum/i))

# print stuff version 2: also include L0 norm of cluster assignment to 000011111 or 111100000, mean of pvalues of D1/D2 genes
diff_000011111 = 0
diff_111100000 = 0
for i in range(0,4):
    if qual[i][1] == '0':
        diff_111100000 += 1
    else:
        diff_000011111 += 1
for i in range(4,9):
    if qual[i][1] == '0':
        diff_000011111 += 1
    else:
        diff_111100000 += 1
sys.stdout.write(str(np.min([diff_000011111,diff_111100000]))+"\t")

sum = 0
i = 1
for entry in qual:
    sum += entry[2]
    i += 1
sys.stdout.write(str(sum/i)+"\t")

sum = 0
i = 1
for entry in qual:
    sum += entry[2]
    i += 1
sys.stdout.write(str(sum/i))
