#!/bin/bash/python

import sys
import os
import numpy as np

remove_these = np.loadtxt(sys.argv[1])
ind = np.zeros(326)
ind[remove_these.astype(int)] = 1
ind = np.where(1-ind)

data_dir = "/data/jessez/Gene_count_datasets/Datasets/batch_all/"
for i in os.listdir(data_dir):
    X = np.loadtxt(data_dir + "/" + i)
    X = np.squeeze(X[ind,:])
    np.savetxt(data_dir+"/"+i+"_NoOutliers",X,delimiter='\t')
