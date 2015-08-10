#!/usr/bin/env python

import sys
import numpy as np
from sklearn import cluster
from clustering_functions import *

def run_clustering(X):
    # Run kmeans clustering. 
    k_means = cluster.KMeans(n_clusters=2)
    k_means.fit(X)
    labels = k_means.labels_
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
