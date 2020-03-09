#!/usr/bin/env python3
#import warnings
#warnings.filterwarnings("ignore")
import umap
#warnings.filterwarnings("always")
from threadpoolctl import threadpool_limits
import numpy as np
from scipy.spatial.distance import squareform
import sys

#random.seed(a=345210)

if __name__=="__main__":
    try:
        variants = np.load(sys.argv[1])
        threads = int(sys.argv[2])
        path = sys.argv[4]

    except IndexError:
            print("Usage <Variant Counts> <Threads> <Path>")
            sys.exit()
    reducer = umap.UMAP(metric='manhattan')
    with threadpool_limits(limits=threads, user_api='blas'):
        variants = variants / np.linalg.norm(variants)
        embedding = reducer.fit_transform(variants)
        np.save(path + "_umap", embedding.astype('float64'), False)
