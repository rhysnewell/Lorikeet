#!/usr/bin/env python3
#import warnings
#warnings.filterwarnings("ignore")
import nimfa
#warnings.filterwarnings("always")
import numpy as np
from scipy.spatial.distance import squareform
import sys

#random.seed(a=345210)

def perform_nmf(array, constraints, k=10, miter=10, estimateRanks='True'):
    # array = [[] for i in filenames]
    # bd = nimfa.Nsnmf(array, seed='nndsvd', rank=k, max_iter=miter, update='euclidean',
    #                 objective='conn')

    bd = nimfa.Bd(array, seed='nndsvd', rank=k, max_iter=miter)

    # bd = nimfa.Pmfcc(array, seed='nndsvd', rank=k, max_iter=miter, theta=constraints)
    bd_fit = bd()
    if estimateRanks == 'True':
        print(bd_fit.fit.rss())

    else:
        print('Rank: %d' % k)
        print('Rss: %5.4f' % bd_fit.fit.rss())
        print('Evar: %5.4f' % bd_fit.fit.evar())
        print('K-L divergence: %5.4f' % bd_fit.distance(metric='kl'))
        print('Sparseness, W: %5.4f, H: %5.4f' % bd_fit.fit.sparseness())
        print('Connectivity', bd_fit.fit.connectivity())

        predictions = bd_fit.fit.predict(prob=True)
        bins = np.array(predictions[0])[0]
        print(predictions)


if __name__=="__main__":
    try:
        pairwise_distances = np.load(sys.argv[4])
        # constraints = np.load(sys.argv[5])
        minRank = int(sys.argv[1])
        estimateRanks = sys.argv[2]
        miter = int(sys.argv[3])
        sample_count = int(sys.argv[6])

    except IndexError:
        print("Usage <Ranks> <Estimate Ranks> <Max Iterations> <Input Pairwise Distance Vector> <Sample Count>")
        sys.exit()

    if sample_count > 0:
        pairwise_distances = squareform(pairwise_distances)
        # constraints = squareform(constraints)

    pairwise_distances = pairwise_distances / np.linalg.norm(pairwise_distances)
    perform_nmf(pairwise_distances, None, minRank, miter, estimateRanks)

