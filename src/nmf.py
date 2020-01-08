import nimfa
import numpy as np
from scipy.spatial.distance import squareform
import sys, os
import ctypes
import ast

def perform_nmf(array, k=10, miter=10, rrange=range(2,20)):
    # array = [[] for i in filenames]
    nsnmf = nimfa.Nsnmf(array, seed='random_vcol', rank=k, max_iter=miter, update='divergence',
                objective='div')
    # lsnmf = nimfa.SepNmf(array, seed='random_vcol', rank=k, n_run=10)
    nsnmf_fit = nsnmf()
    best_rank = nsnmf.estimate_rank(rank_range=rrange)
    for rank, values in best_rank.items():
        print('Rank: %d' % rank)
        print('Rss: %5.4f' % values['rss'])
        print('Evar: %5.4f' % values['evar'])
        print('K-L: %5.4f' % values['kl'])
    # print(best_rank)
    print('Rss: %5.4f' % nsnmf_fit.fit.rss())
    print('Evar: %5.4f' % nsnmf_fit.fit.evar())
    print('K-L divergence: %5.4f' % nsnmf_fit.distance(metric='kl'))
    print('Sparseness, W: %5.4f, H: %5.4f' % nsnmf_fit.fit.sparseness())

    # predictions = nsnmf_fit.fit.predict(prob=True)
    # bins = np.array(predictions[0])[0]
    # bin_dict = {}
    # prob_idx = 0
    # for (bin_id, contig_name) in zip(bins, col_ids):
    #     if predictions[1][prob_idx] >= 0.5:
    #         try:
    #             bin_dict[bin_id].append(contig_name)
    #         except KeyError:
    #             bin_dict[bin_id] = [contig_name]

    # return nsnmf_fit, bin_dict, best_rank

if __name__=="__main__":
    try:
        pairwise_distances = sys.argv[4]
        minrank = int(sys.argv[1])
        maxrank = int(sys.argv[2])
        miter = int(sys.argv[3])
    except IndexError:
        print("Usage <Ranks> <Max Iterations> <Rank Range> <Input Pairwise Distance Vector>")
        sys.exit()

    condensed_vec = []
    with open(pairwise_distances) as f:
        for line in f:
            condensed_vec = ast.literal_eval(line)
            #condensed_vec = line
    condensed_vec = np.array(condensed_vec)
    square_vec = squareform(condensed_vec)

    perform_nmf(square_vec, minrank, miter, range(minrank, maxrank+1))

