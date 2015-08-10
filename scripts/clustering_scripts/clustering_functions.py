import numpy as np
import scipy as sp

def cluster_and_account_for_outliers(X,n,run_clustering,loc,cells):
    labels = run_clustering(X)
    n_0 = np.sum(labels == 0)
    n_1 = np.sum(labels == 1)
    ind = np.array(range(0,len(X)))
    # if there are 3 or less samples in a cluster, they're probably outliers!                                                                                                                             
    while np.min([n_0,n_1]) <= 3 and n_0 + n_1 >= n/2:
        # track outlier indices                                                                                                                                                                           
        if n_0 < n_1:
            outlier_ind = np.delete(range(0,len(labels)),np.where(labels == 1))
        else:
            outlier_ind = np.delete(range(0,len(labels)),np.where(labels == 0))
        # remove outliers and repeat clustering                                                                                                                                                           
        X = np.delete(X,outlier_ind,axis=0)
        ind = np.delete(ind,outlier_ind)
        labels = run_clustering(X)
        n_0 = np.sum(labels == 0)
        n_1 = np.sum(labels == 1)

    # save outlier cells (if any exist)
    if len(ind) < int(n):
        save_outliers(ind,loc,cells)

    # generate random labels for the outliers (keep same proportions as generated labels)
    p_0 = np.sum(labels == 0)/float(len(labels))
    final_labels = np.zeros(n)
    ii = 0
    for i in range(0,int(n)):
        if i in ind:
            final_labels[i] = labels[ii]
            ii += 1
        else:
            if np.random.uniform() < p_0:
                final_labels[i] = 0
            else:
                final_labels[i] = 1

    return final_labels.astype(int)

def save_outliers(ind,loc,cells):
    outlier_cells = np.delete(cells,ind).astype(int)
    np.savetxt(loc,outlier_cells,fmt='%s')
