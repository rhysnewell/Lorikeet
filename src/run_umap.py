#!/usr/bin/env python3
#import warnings
#warnings.filterwarnings("ignore")
import umap
#warnings.filterwarnings("always")
from threadpoolctl import threadpool_limits
import numpy as np
import matplotlib.pyplot as plt
import sys

#random.seed(a=345210)

if __name__=="__main__":
    try:
        variants = np.load(sys.argv[1])
        threads = int(sys.argv[2])

    except IndexError:
            print("Usage <Variant Counts> <Threads>")
            sys.exit()
    reducer = umap.UMAP(metric='manhattan')
    with threadpool_limits(limits=threads, user_api='blas'):
        variants = variants / np.linalg.norm(variants)
        embedding = reducer.fit_transform(variants)
        plt.scatter(embedding[:, 0], embedding[:, 1])
        plt.gca().set_aspect('equal', 'datalim')
        plt.title('UMAP projection of variants', fontsize=24)
        plt.savefig("umap_projection.png")

        np.save(sys.argv[1] + "-umap", embedding.astype('float64'), False)
