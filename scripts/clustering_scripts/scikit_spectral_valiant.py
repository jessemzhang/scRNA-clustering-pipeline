#!/usr/bin/env python

import sys
import numpy as np
from sklearn import cluster
import sklearn.metrics
from clustering_functions import *

def valiant_dist(x,y):
    x = x.astype(float)
    y = y.astype(float)
    a = x + y
    b = x - y
    c = (a != 0)
    # for entries that have x + y = 0, the contribution should be -1)
    D = (np.sum((np.square(b[c])-a[c])/a[c])-(len(x)-np.sum(c)))
    return D

def compute_affinity_matrix(X):
    # Find the L2 norm between each two pair of points..
    m = X.shape[0]
    D = np.zeros([m,m])
    for i in range(0,m):
        for j in range(0,m):
            D[i,j] = valiant_dist(X[i,:],X[j,:])
    # Apply gaussian kernel to transform distance matrix into similarity matrix
    delta = 1e6
    return np.exp(- D ** 2 / (2. * delta ** 2))

def run_clustering(X):
    # Run spectral clustering. 
    D = compute_affinity_matrix(X)
    spectral = cluster.SpectralClustering(n_clusters=2,affinity='precomputed')
    spectral.fit(D)
    labels = spectral.labels_
    return labels

if __name__ == "__main__":

    # LOAD DATA                                                                                  
    X = np.loadtxt(sys.argv[1])
    loc = sys.argv[2]
    cells = np.loadtxt(sys.argv[3])
    
    n = float(np.shape(X)[0]) # number of samples                                                       

    labels = cluster_and_account_for_outliers(X,n,run_clustering,loc,cells)
    for label in labels:
        print label

