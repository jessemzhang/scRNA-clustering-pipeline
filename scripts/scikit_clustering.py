#!/usr/bin/env python
#compatibility with python 3

import sys
from time import time
import numpy as np;
import matplotlib.pyplot as plt;
from pylab import *
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble, lda,
                     random_projection)


# LOAD DATA
X_data = np.loadtxt(sys.argv[1])
#X_data = np.transpose(X_data)
print(shape(X_data))


# FIT MODEL
model = manifold.TSNE(n_components=2, random_state=0);
print('fitting...')
X_fit = model.fit_transform(X_data);

# SAVE LABELS
print('saving...')
np.savetxt(sys.argv[2],X_fit)

