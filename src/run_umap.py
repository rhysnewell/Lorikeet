#!/usr/bin/env python3
#import warnings
#warnings.filterwarnings("ignore")
import umap
#warnings.filterwarnings("always")
from threadpoolctl import threadpool_limits
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import statistics
import numba
from sklearn.decomposition import PCA

#random.seed(a=345210)

"""
This Script is a placeholder for running umap from Rust using the python module, umap-learn
Eventually we will migrate umap over to rust
"""

@numba.njit()
def coda_rho(a, b):
    var_a = statistics.variance(a)
    var_b = statistics.variance(b)
    covar = np.cov(a, b)[0][1]

    rho = 2*covar / (var_a + var_b)

    # rho *= -1
    # rho += 1
    return covar

if __name__=="__main__":
    try:
        variants = np.load(sys.argv[1])
        threads = int(sys.argv[2])
        min_dist = float(sys.argv[3])
        spread = float(sys.argv[4])
        n_neighbours = int(sys.argv[5])
        metric = sys.argv[6]

    except IndexError:
        print("Usage <Variant Counts> <Threads> <Min Dist> <Spread> <No. Neighbours> <Metric>")
        sys.exit()
    reducer = umap.UMAP(metric=metric, min_dist=min_dist, spread=spread, n_neighbors=n_neighbours, n_components=2)
    with threadpool_limits(limits=threads, user_api='blas'):
        variants = variants / np.linalg.norm(variants, ord=2, axis=1, keepdims=True)
        embedding = reducer.fit_transform(variants)
        plt.scatter(embedding[:, 0], embedding[:, 1], alpha=0.33)
        plt.gca().set_aspect('equal', 'datalim')
        plt.title('UMAP projection of variants', fontsize=24)
        plt.savefig("umap_projection.png")

        np.save(sys.argv[1] + "-umap", embedding.astype('float64'), False)
