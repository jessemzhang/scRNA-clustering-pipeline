#!/usr/bin/env python

# This script requires 3 input arguments (in this order): the data matrix, the labels, and gene names. It uses random forests to determine which genes are the most important for classifying the data into the two labels. Note that if the data matrix is m-by-n, the label vector is length m, and the gene names vector is length n. (m observations). min(100,round(m/2)) gene names are reported.

from sklearn.ensemble import RandomForestClassifier;
import numpy as np
import scipy as sp
import sys

def run(X,Y):
    clf = RandomForestClassifier();
    clf.fit(X,Y);
    return clf.feature_importances_;

def select_genes_using_RF(X,Y):
    # find most important features using random forest
    n = 100; # number of times to repeat random forest
    m = int(np.round(X.shape[1]/2)); # number of features to output at the end
    if m > n:
        m = n
    
    rankings = np.zeros((X.shape[1],1))
    for i in range(0,n):    
        importances = run(X,Y)
        most_important = np.flipud(np.argsort(importances))
        rankings[most_important[0:m]] += 1
    
    most_important = np.flipud(np.transpose(np.argsort(np.transpose(rankings))))[0:m]
    return most_important

# load data
X = np.loadtxt(sys.argv[1]);
Y = np.loadtxt(sys.argv[2],dtype=str)
gene_names = np.loadtxt(sys.argv[3],dtype=str)

important_genes = gene_names[select_genes_using_RF(X,Y)]

# f = open('important_genes.txt','w')
for gene in important_genes:
    # print >> f, gene[0]
    print gene[0]
