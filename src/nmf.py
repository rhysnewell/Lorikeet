#!/usr/bin/env python3
#import warnings
#warnings.filterwarnings("ignore")
import nimfa
#warnings.filterwarnings("always")
import numpy as np
from scipy.spatial.distance import squareform
import sys

#random.seed(a=345210)

def perform_nmf(array, k=10, miter=10, estimateRanks='True'):
    # array = [[] for i in filenames]
    bd = nimfa.Nsnmf(array, seed='nndsvd', rank=k, max_iter=miter, update='euclidean',
                    objective='conn')
    bd_fit = bd()
    if estimateRanks == 'True':
        print(bd_fit.fit.rss())

    else:

        print('Rss: %5.4f' % bd_fit.fit.rss())
        print('Evar: %5.4f' % bd_fit.fit.evar())
        print('K-L divergence: %5.4f' % bd_fit.distance(metric='kl'))
        print('Sparseness, W: %5.4f, H: %5.4f' % bd_fit.fit.sparseness())
        print('Connectivity', bd_fit.fit.connectivity())

        predictions = bd_fit.fit.predict(prob=True)
        bins = np.array(predictions[0])[0]
        print(predictions)
    #print('Evar: %5.4f' % bd_fit.fit.evar())
    #print('K-L divergence: %5.4f' % bd_fit.distance(metric='kl'))
    #print('Sparseness, W: %5.4f, H: %5.4f' % bd_fit.fit.sparseness())

    # lsnmf = nimfa.SepNmf(array, seed='random_vcol', rank=k, n_run=10)
    #estimated_ranks = bd.estimate_rank(rank_range=rrange, n_run=miter//2)

    #best_rank = 0
    #best_rss = None
    ## Choose the best rank by finding the first inflection point in the RSS values
    ## Quick method is to find first point when RSS is at lowest value
    #for rank, values in estimated_ranks.items():
    #    print('Rank: %d' % rank)
    #    print('Rss: %5.4f' % values['rss'])
    #    print('Evar: %5.4f' % values['evar'])
    #    print('K-L: %5.4f' % values['kl'])
    #    if best_rss is None or values['rss'] < best_rss:
    #        best_rank = rank
    #        best_rss = values['rss']
    #    elif values['rss'] >= best_rss:
    #        break

    #print("%d strains estimated" % best_rank)
    #bd = nimfa.Nmf(array, seed='nndsvd', rank=best_rank, max_iter=miter, update='euclidean',
    #                                     objective='fro')

    #bd_fit = bd()
    #print('Rss: %5.4f' % bd_fit.fit.rss())
    #print('Evar: %5.4f' % bd_fit.fit.evar())
    #print('K-L divergence: %5.4f' % bd_fit.distance(metric='kl'))
    #print('Sparseness, W: %5.4f, H: %5.4f' % bd_fit.fit.sparseness())



    #predictions = bd_fit.fit.predict(prob=True)
    #bins = np.array(predictions[0])[0]
    #print(predictions)

    # return nsnmf_fit, bin_dict, best_rank

if __name__=="__main__":
    try:
        pairwise_distances = np.load(sys.argv[4])
        minRank = int(sys.argv[1])
        estimateRanks = sys.argv[2]
        miter = int(sys.argv[3])
    except IndexError:
        print("Usage <Ranks> <Estimate Ranks> <Max Iterations> <Input Pairwise Distance Vector>")
        sys.exit()

    #condensed_vec = []
    #with open(pairwise_distances) as f:
    #    for line in f:
    #        condensed_vec = ast.literal_eval(line)
    #        #condensed_vec = line
    #condensed_vec = np.array(condensed_vec)
    #pairwise_distances = squareform(pairwise_distances)

    perform_nmf(pairwise_distances, minRank, miter, estimateRanks)

