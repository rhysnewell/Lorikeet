#!/usr/bin/env python3
#import warnings
#warnings.filterwarnings("ignore")
import nimfa
#warnings.filterwarnings("always")
from threadpoolctl import threadpool_limits
import numpy as np
from scipy.spatial.distance import squareform
import sys

#random.seed(a=345210)

def perform_nmf(array, constraints, k=10, miter=10, estimateRanks='True', path='/tmp/lorikeet-pred'):
    # array = [[] for i in filenames]
    # mf = nimfa.Nsnmf(array, seed='nndsvd', rank=k, max_iter=miter, update='euclidean',
    #                 objective='conn')
    if constraints is None:
        converged = False
        count = 0
        while converged is False or count >= 10:
            try:
                mf = nimfa.Nmf(array, seed='nndsvd', rank=k, max_iter=miter, update='euclidean',
                                objective='conn')
                # mf = nimfa.Psmf(array, seed=None, rank=k, max_iter=miter)
                mf_fit = mf()
                converged  = True
            except LinAlgError:
                count += 1

        if estimateRanks == 'True':
            print(mf_fit.fit.rss())

        else:
            if count >= 10:
                converged = False
                while converged is False or k >= 1:
                    try:
                        mf = nimfa.Nmf(array, seed='nndsvd', rank=k, max_iter=miter, update='euclidean',
                                        objective='conn')
                        # mf = nimfa.Psmf(array, seed=None, rank=k, max_iter=miter)
                        mf_fit = mf()
                        converged  = True
                    except LinAlgError:
                        k -= 1
            else if converged:
                print('Rank: %d' % k)
                print('Rss: %5.4f' % mf_fit.fit.rss())
                print('Evar: %5.4f' % mf_fit.fit.evar())
                print('K-L divergence: %5.4f' % mf_fit.distance(metric='kl'))
                print('Sparseness, W: %5.4f, H: %5.4f' % mf_fit.fit.sparseness())
                print('Connectivity', mf_fit.fit.connectivity())
                print('Score', mf_fit.fit.select_features())

                predictions = mf_fit.fit.predict(prob=True)
                bins = np.array(predictions[0])[0]
                # bins = np.array(predictions[0])[0].reshape((len(bins), 1))
                new_pred = np.column_stack((bins,
                                            np.array(predictions[1]),
                                            np.array(mf_fit.fit.select_features())))
                print(path)

                np.save(path, new_pred.astype('float32'), False)
            else:
                print("Failed to converge SVD, Try again with stricter variant calling")
    else:
        mf = nimfa.Pmfcc(array, seed='nndsvd', rank=k, max_iter=miter, theta=constraints)
        mf_fit = mf()
        if estimateRanks == 'True':
            print(mf_fit.distance())



if __name__=="__main__":
    try:
        pairwise_distances = np.load(sys.argv[4])
        minRank = int(sys.argv[1])
        estimateRanks = sys.argv[2]
        miter = int(sys.argv[3])
        sample_count = int(sys.argv[6])
        threads = int(sys.argv[7])
        try:
            # constraints = np.load(sys.argv[5])
            constraints = None
        except FileNotFoundError:
            constraints = None

    except IndexError:
            print("Usage <Ranks> <Estimate Ranks> <Max Iterations> <Input Pairwise Distance Vector> <BLANK> <Sample Count> <Threads>")
            sys.exit()

    with threadpool_limits(limits=threads, user_api='blas'):
        if sample_count > 0:
            pairwise_distances = squareform(pairwise_distances)
            if constraints is not None:
                constraints = squareform(constraints)

        pairwise_distances = pairwise_distances / np.linalg.norm(pairwise_distances)
        perform_nmf(pairwise_distances, constraints, minRank, miter, estimateRanks, sys.argv[4])

