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

#random.seed(a=345210)

@numba.njit()
def coda_rho(a, b):
    var_a = statistics.variance(a)
    var_b = statistics.variance(b)
    covar = np.cov(a, b)[0][1]

    rho = 2*covar / (var_a + var_b)

    rho *= -1
    rho += 1
    return rho

if __name__=="__main__":
    try:
        variants = np.load(sys.argv[1])
        threads = int(sys.argv[2])

    except IndexError:
        print("Usage <Variant Counts> <Threads>")
        sys.exit()
    reducer = umap.UMAP(metric="mahalanobis", min_dist=0.15, spread=0.5, n_neighbors=int(0.1*np.size(variants, 0)))
    with threadpool_limits(limits=threads, user_api='blas'):
        variants = variants / np.linalg.norm(variants, ord=2, axis=1, keepdims=True)
        embedding = reducer.fit_transform(variants)
        plt.scatter(embedding[:, 0], embedding[:, 1], alpha=0.33)
        plt.gca().set_aspect('equal', 'datalim')
        plt.title('UMAP projection of variants', fontsize=24)
        plt.savefig("umap_projection.png")

        np.save(sys.argv[1] + "-umap", embedding.astype('float64'), False)
