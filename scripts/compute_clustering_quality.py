# Given a data matrix and labels, generate a measure of clustering quality using some distance metric

import sys
import numpy as np
import sklearn.metrics

x = np.loadtxt(sys.argv[1])
y = np.loadtxt(sys.argv[2],dtype=str)
n = len(x)

d = sklearn.metrics.pairwise.euclidean_distances(x)
# d = sklearn.metrics.pairwise.cosine_similarity(x)

# labels
labels = np.unique(y)

# average distance between points within clusters (assuming 0 along diagonals)
d_0 = d[y == labels[0],]
d_0 = d_0[:,y == labels[0]]
d_1 = d[y == labels[1],]
d_1 = d_1[:,y == labels[1]]
d_within = (np.sum(d_0)+np.sum(d_1))/float(len(d_0)**2-len(d_0)+len(d_1)**2-len(d_1))

# average distance between points in different clusters
d_0_1 = d[y == labels[0],]
d_0_1 = d_0_1[:,y == labels[1]]
d_between = np.sum(d_0_1)/(np.sum(y == labels[0])*np.sum(y == labels[1]))

# print!
print d_between/d_within
