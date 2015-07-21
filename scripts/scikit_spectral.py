#!/usr/bin/env python

import sys
import numpy as np
from sklearn import cluster
import sklearn.metrics

def compute_affinity_matrix(X):
    # valid metrics are: 'rbf', 'sigmoid', 'polynomial', 'poly', 'linear', 'cosine'
    D = sklearn.metrics.pairwise.pairwise_kernels(X,metric='rbf')
    return D

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
    labels = run_clustering(X)
    for label in labels:
        print label

