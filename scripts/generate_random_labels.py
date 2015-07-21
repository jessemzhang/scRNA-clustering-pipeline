#!/usr/bin/python

import sys
import numpy as np

#num_entries = int(sys.stdin.read())

num_entries = int(sys.argv[1])
ratio = int(sys.argv[2])
x = np.random.uniform(0,1,(num_entries,1))

f = open('../labels/labels_random.txt','w')
i = 0
for entry in x:
    if x[i] > ratio:
        print >> f, '1'
    else:
        print >> f, '0'
    i += 1
