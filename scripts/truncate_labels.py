#!/bin/bash/python

import numpy as np
import sys

cells=np.loadtxt(sys.argv[2],dtype=str)
labels=np.loadtxt(sys.argv[1])

# dictionary for mapping all cells         
all_cells = np.loadtxt('/data/jessez/Gene_count_datasets/Cells/cells_all.txt',dtype=str);
dictionary = dict(zip(all_cells, range(0,len(all_cells))))

# cut down indices to only those in the data matrix
ind = np.array(map(lambda x:dictionary[x],cells))
labels_trunc = labels[ind]

# print the labels..
for label in labels_trunc:
    print int(label)
