#!/usr/bin/env python

import sys
import numpy as np
from sklearn import cluster

def run_clustering(X):
    # Run kmeans clustering. 

    k_means = cluster.KMeans(n_clusters=2)
    k_means.fit(X)
    labels = k_means.labels_
    return labels

if __name__ == "__main__":

    # LOAD DATA
    X = np.loadtxt(sys.argv[1])

    labels = run_clustering(X)
    for label in labels:
        print label

