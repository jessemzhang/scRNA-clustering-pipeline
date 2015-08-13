#!/usr/bin/env python

# This script computes the number of disagreeing labels between truth labels and another set of labels. The only input is:
# -labels 
# -cells

import sys
import numpy as np

# dictionary for mapping all cells
all_cells = np.loadtxt('/data/jessez/Gene_count_datasets/Cells/cells_all.txt',dtype=str);
dictionary = dict(zip(all_cells, range(0,len(all_cells))))

labels = np.loadtxt(sys.argv[1])
labels_d1 = np.loadtxt('../labels/labels_d1_truth.txt')
labels_d2 = np.loadtxt('../labels/labels_d2_truth.txt')

cells = np.loadtxt(sys.argv[2],dtype=str)

# cut down the truth indices to only those of cells in labels
ind = np.array(map(lambda x:dictionary[x],cells))
labels_d1 = labels_d1[ind]
labels_d2 = labels_d2[ind]

d1_check = labels[np.where(labels_d1)] # should all be 1 or 0 if good labels
d2_check = labels[np.where(labels_d2)] # should all be 1 or 0 if good labels

# compute hamming distance..

# if 0 is d1 and 1 is d2, then # errors is:
err_0_is_d1 = np.sum(d1_check == 1)+np.sum(d2_check == 0)
# otherwise, 1 is d1 and 2 is d2, and the # errors is:
err_1_is_d1 = np.sum(d1_check == 0)+np.sum(d2_check == 1)

# output the error rate (choose min)
print np.min([err_0_is_d1,err_1_is_d1])/float(np.shape(d1_check)[0]+np.shape(d2_check)[0])
