#!/usr/bin/env python

import sys
import numpy as np
from sklearn import cluster

def run_clustering(X):
    # Run kmeans clustering. To help this function generate two roughly equal-sized clusters, have k-means create 3 clusters. hopefully the two smaller clusters together have the same amount of samples as the bigger cluster.. 

    k_means = cluster.KMeans(n_clusters=3)
    k_means.fit(X)
    labels = k_means.labels_

    a = np.unique(labels,return_counts=True)
    labels[labels==np.argmax(a[1])] = -1
    labels[labels!=-1] = 0
    labels[labels==-1] = 1
    return labels

if __name__ == "__main__":

    # LOAD DATA
    X = np.loadtxt(sys.argv[1])

    labels = run_clustering(X)
    for label in labels:
        print label

