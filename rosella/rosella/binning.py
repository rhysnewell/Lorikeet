#!/usr/bin/env python
###############################################################################
# binning.py - A fast binning algorithm spinning off of the methodology of
#              Lorikeet
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################

__author__ = "Rhys Newell"
__copyright__ = "Copyright 2020"
__credits__ = ["Rhys Newell"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Rhys Newell"
__email__ = "rhys.newell near hdr.qut.edu.au"
__status__ = "Development"

###############################################################################
# System imports
import sys
import argparse
import logging
import os
import shutil
import datetime

# Function imports
import numpy as np
from numba import njit
import multiprocessing as mp
import pandas as pd
import hdbscan
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from Bio import SeqIO
import skbio.stats.composition
import umap
from itertools import permutations

# self imports
import rosella.metrics as metrics

# Set plotting style
sns.set(style='white', context='notebook', rc={'figure.figsize': (14, 10)})

# Debug
debug = {
    1: logging.CRITICAL,
    2: logging.ERROR,
    3: logging.WARNING,
    4: logging.INFO,
    5: logging.DEBUG
}

###############################################################################
############################### - Exceptions - ################################


class BadTreeFileException(Exception):
    pass


###############################################################################                                                                                                                      [44/1010]
################################ - Functions - ################################


def phelp():
    print("""
Usage:
rosella.py [SUBCOMMAND] ..

Subcommands:
fit
""")


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def write_contig(contig, assembly, f):
    seq = assembly[contig]
    fasta = ">" + seq.id + '\n'
    fasta += str(seq.seq) + '\n'
    f.write(fasta)

@njit
def index(array, item):
    for idx, val in np.ndenumerate(array):
        if val == item:
            return idx

###############################################################################
################################ - Classes - ##################################


class CustomHelpFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if '%(default)' not in action.help:
            if action.default != '' and \
               action.default != [] and \
               action.default != None \
               and action.default != False:
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [
                        argparse.OPTIONAL, argparse.ZERO_OR_MORE
                    ]

                    if action.option_strings or action.nargs in defaulting_nargs:

                        if '\n' in h:
                            lines = h.splitlines()
                            lines[0] += ' (default: %(default)s)'
                            h = '\n'.join(lines)
                        else:
                            h += ' (default: %(default)s)'
        return h

    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])


class Binner():
    def __init__(
            self,
            count_path,
            output_prefix,
            assembly,
            scaler="clr",
            n_neighbors=20,
            min_dist=0.1,
            n_components=2,
            random_state=42,
            min_cluster_size=100,
            min_contig_size=2500,
            min_samples=50,
            prediction_data=True,
            cluster_selection_method="eom",
            precomputed=False,
            hdbscan_metric="euclidean",
            metric = 'aggregate',
            threads=8,
    ):
        self.pool = mp.Pool(threads)
        # Open up assembly
        self.assembly = SeqIO.to_dict(SeqIO.parse(assembly, "fasta"))

        ## Set up clusterer and UMAP
        self.path = output_prefix
        self.coverage_table = pd.read_csv(count_path, sep='\t')
        self.large_contigs = self.coverage_table[self.coverage_table["contigLen"] >= min_contig_size]
        # self.small_contigs = self.coverage_table[(1000 <= self.coverage_table["contigLen"]) & (self.coverage_table["contigLen"] < min_contig_size)]

        # If there are enough contigs of that size
        if self.large_contigs.shape[0] > 100:
            self.depths = self.large_contigs.iloc[:,3:]
            # self.small_depths = self.small_contigs.iloc[:,3:]
        else: # Otherwise we'll just use a smaller value
            self.large_contigs = self.coverage_table[self.coverage_table["contigLen"] >= 1000]
            # self.small_contigs = self.coverage_table[(500 < self.coverage_table["contigLen"]) & (self.coverage_table["contig_len"] < 1000)]
            self.depths = self.large_contigs.iloc[:,3:]
            # self.small_depths = self.small_contigs.iloc[:,3:]

        # if self.depths.shape[1] > 2:
        self.depths = self.depths[self.depths.columns[::2]]
            # self.small_depths = self.small_depths[self.small_depths.columns[::2]]

        logging.info("Calculating TNF values")
        ## Add the TNF values
        tnf_dict = {}
        for (idx, contig) in enumerate(self.coverage_table.iloc[:,0]):
            seq = self.assembly[contig].seq
            for s in [str(seq).upper(),
                      str(seq.reverse_complement()).upper()]:
                # For di, tri and tetranucleotide counts, we loop over the
                # sequence and its reverse complement, until we're near the end:
                for i in range(len(s[:-4])):
                    tetra = s[i:i + 4]
                    try:
                        tnf_dict[str(tetra)][idx] += 1
                    except:
                        tnf_dict[str(tetra)] = [0] * self.coverage_table.iloc[:,0].values.shape()[0]
                        tnf_dict[str(tetra)][idx] += 1

        for (tnf, vector) in tnf_dict.items():
            vector_sum = sum(vector)
            freqs = [x/vector_sum for x in vector]
            self.depths[tnf] = freqs

        ## Scale the data
        if scaler.lower() == "minmax":
            self.depths = MinMaxScaler().fit_transform(self.depths)
            self.small_depths = MinMaxScaler().fit_transform(self.small_depths)
        elif scaler.lower() == "clr":
            self.depths = skbio.stats.composition.clr(self.depths + 1)
            # self.small_depths = skbio.stats.composition.clr(self.small_depths + 1)
        elif scaler.lower() == "none":
            pass

        if n_neighbors >= int(self.depths.shape[0] * 0.5):
            n_neighbors = max(int(self.depths.shape[0] * 0.5), 2)

        if n_components > self.depths.shape[1]:
            n_components = self.depths.shape[1]

        if metric in ['aggregate', 'rho', 'phi', 'phi_dist']:
            self.reducer = umap.UMAP(
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                n_components=n_components,
                random_state=random_state,
                spread=1,
                metric=getattr(metrics, metric)
            )
        else:
            self.reducer = umap.UMAP(
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                n_components=n_components,
                random_state=random_state,
                spread=1,
                metric=metric
            )


        if min_cluster_size > self.depths.shape[0] * 0.1:
            min_cluster_size = max(int(self.depths.shape[0] * 0.1), 2)
            min_samples = max(int(min_cluster_size * 0.1), 2)

        if precomputed:
            metric = "precomputed"
            prediction_data = False

        self.clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            # min_samples=min_samples,
            prediction_data=prediction_data,
            cluster_selection_method=cluster_selection_method,
            metric=hdbscan_metric,
        )

    def fit_transform(self):
        ## Calculate the UMAP embeddings
        logging.info("Running UMAP - %s" % self.reducer)
        self.embeddings = self.reducer.fit_transform(self.depths)
        # self.small_embeddings = self.reducer.transform(self.small_depths)

    def cluster(self):
        ## Cluster on the UMAP embeddings and return soft clusters
        try:
            logging.info("Running HDBSCAN - %s" % self.clusterer)
            self.clusterer.fit(self.embeddings)
            self.soft_clusters = hdbscan.all_points_membership_vectors(
                self.clusterer)
            self.soft_clusters_capped = np.array([np.argmax(x) for x in self.soft_clusters])
            # self.small_labels, self.small_strengths = hdbscan.approximate_predict(self.clusterer, self.small_embeddings)
        except:
            ## Likely integer overflow in HDBSCAN
            ## Try reduce min samples
            self.clusterer = hdbscan.HDBSCAN(
                min_cluster_size=max(int(self.depths.shape[0] * 0.01), 2),
                min_samples=max(int(self.depths.shape[0] * 0.005), 2),
                prediction_data=True,
                cluster_selection_method="eom",
            )
            logging.info("Retrying HDBSCAN - %s" % self.clusterer)
            self.clusterer.fit(self.embeddings)
            self.soft_clusters = hdbscan.all_points_membership_vectors(
                self.clusterer)
            self.soft_clusters_capped = np.array([np.argmax(x) for x in self.soft_clusters])
            # self.small_labels, self.small_strengths = hdbscan.approximate_predict(self.clusterer, self.small_embeddings)


    def cluster_distances(self):
        ## Cluster on the UMAP embeddings and return soft clusters
        try:
            logging.info("Running HDBSCAN - %s" % self.clusterer)
            self.clusterer.fit(self.depths)
        except:
            ## Likely integer overflow in HDBSCAN
            ## Try reduce min samples
            self.clusterer = hdbscan.HDBSCAN(
                min_cluster_size=max(int(self.depths.shape[0] * 0.01), 2),
                min_samples=max(int(self.depths.shape[0] * 0.005), 2),
                prediction_data=True,
                cluster_selection_method="precomputed",
            )
            logging.info("Retrying HDBSCAN - %s" % self.clusterer)

            self.clusterer.fit(self.depths)

    def plot(self):
        logging.info("Generating UMAP plot with labels")

        label_set = set(self.clusterer.labels_)
        color_palette = sns.color_palette('Paired', len(label_set))
        cluster_colors = [
            color_palette[x] if x >= 0 else (0.5, 0.5, 0.5) for x in self.soft_clusters_capped
        ]

        # small_cluster_colors = [
        #     color_palette[x] if x >= 0 else (0.5, 0.5, 0.5) for x in self.small_labels
        # ]

        cluster_member_colors = [
            sns.desaturate(x, p) for x, p in zip(cluster_colors, self.clusterer.probabilities_)
        ]
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ## Plot large contig membership
        ax.scatter(self.embeddings[:, 0],
                   self.embeddings[:, 1],
                   s=7,
                   linewidth=0,
                   c=cluster_member_colors,
                   alpha=0.7)

        # ## Plot small contig membership
        # ax.scatter(self.small_embeddings[:, 0],
        #            self.small_embeddings[:, 1],
        #            s=5,
        #            linewidth=0,
        #            c=small_cluster_colors,
        #            alpha=0.7)
        # ax.add_artist(legend)
        plt.gca().set_aspect('equal', 'datalim')
        plt.title('UMAP projection of contigs', fontsize=24)
        plt.savefig(self.path + '/UMAP_projection_with_clusters.png')

    def plot_distances(self):
        label_set = set(self.clusterer.labels_)
        self.clusterer.condensed_tree_.plot(
            select_clusters=True,
            selection_palette=sns.color_palette('deep', len(label_set)),
        )
        plt.title('Hierarchical tree of clusters', fontsize=24)
        plt.savefig(self.path + '/cluster_hierarchy.png')

    def labels(self):
        try:
            return self.soft_clusters.astype('int8')
        except AttributeError:
            return self.clusterer.labels_.astype('int8')

    def bin_contigs(self, assembly_file, min_bin_size=200000):
        logging.info("Binning contigs...")

        # initialize bin dictionary Label: Vec<Contig>
        self.bins = {}
        for (idx, label) in enumerate(self.clusterer.labels_):
            if label != -1:
                try:
                    self.bins[label].append(idx)
                except KeyError:
                    self.bins[label] = [idx]
            # elif len(self.assembly[self.coverage_table.iloc[idx, 0]].seq) >= min_bin_size:
            #     try:
            #         self.bins[label].append(idx)
            #     except KeyError:
            #         self.bins[label] = [idx]
            else:
                soft_label = self.soft_clusters_capped[idx]
                try:
                    self.bins[soft_label].append(idx)
                except KeyError:
                    self.bins[soft_label] = [idx]

        # ## Bin out small contigs
        # for (idx, label) in enumerate(self.small_labels):
        #     if label != -1:
        #         try:
        #             self.bins[label].append(idx)
        #         except KeyError:
        #             self.bins[label] = [idx]

        logging.info("Merging bins...")
        for bin in list(self.bins):
            if bin != -1:
                self.pool.apply_async(self.spawn_merge, args=(bin, min_bin_size))

        self.pool.close()
        self.pool.join()  # postpones the execution of next line of code until all processes in the queue are done.

        logging.info("Writing bins...")
        for (bin, contigs) in self.bins.items():
            if bin != -1:
                # Calculate total bin size and check if it is larger than min_bin_size
                bin_length = sum([len(self.assembly[self.coverage_table.iloc[idx, 0]].seq) for idx in contigs])
                if bin_length >= min_bin_size:
                    with open(self.path + '/rosella_bin.' + str(bin) + '.fna', 'w') as f:
                        for idx in contigs:
                            contig = self.coverage_table.iloc[idx, 0]
                            write_contig(contig, self.assembly, f)

            else:
                # Get final bin value
                bin_max = max(self.bins.keys()) + 1
                # Rescue any large unbinned contigs and put them in their own cluster
                for idx in contigs:
                    contig = self.coverage_table.iloc[idx, 0]
                    if len(self.assembly[self.coverage_table.iloc[idx, 0]].seq) >= min_bin_size:
                        with open(self.path + '/rosella_bin.' + str(bin_max) + '.fna', 'w') as f:
                            write_contig(contig, self.assembly, f)
                        bin_max += 1

    def spawn_merge(self, bin, min_bin_size=300000):
        # Calculate total bin size and check if it is larger than min_bin_size
        contigs = self.bins[bin]
        bin_length = sum([len(self.assembly[self.coverage_table.iloc[idx, 0]].seq) for idx in contigs])
        if bin_length < min_bin_size:
            # reassign these contigs
            for idx in contigs:
                self.pool.apply_async(self.merge_bins, args=(idx), callback=self.collect_merge)

    def merge_bins(self, idx):
        soft_clusters = self.soft_clusters[idx]
        second_max = sorted(soft_clusters, reverse=True)[1]
        try:
            next_label = index(soft_clusters, second_max)[0]
        except IndexError:
            next_label = -1

        return (next_label, idx)


    def collect_merge(self, result):
        try:
            self.bins[result[0]].append(result[1])
            # self.bins[bin].remove(idx)
        except KeyError:
            self.bins[result[0]] = [result[1]]


if __name__ == '__main__':

    ############################ ~ Main Parser ~ ##############################
    main_parser = argparse.ArgumentParser(prog='cluster',
                                          formatter_class=CustomHelpFormatter,
                                          add_help=False)
    main_parser.add_argument('--version',
                             action='version',
                             version=__version__,
                             help='Show version information.')
    main_parser.add_argument(
        '--verbosity',
        help=
        '1 = critical, 2 = error, 3 = warning, 4 = info, 5 = debug. Default = 4 (logging)',
        type=int,
        default=4)
    main_parser.add_argument('--log',
                             help='Output logging information to file',
                             default=False)
    subparsers = main_parser.add_subparsers(help="--", dest='subparser_name')

    ########################## ~ sub-parser ~ ###########################
    input_options = subparsers.add_parser(
        'fit',
        description='Perform UMAP and then HDBSCAN on array of variant depths',
        formatter_class=CustomHelpFormatter,
        epilog='''
                            ~ fit ~
How to use fit:

rosella.py fit --input coverm_output.tsv --assembly scaffolds.fasta

''')
    ## Main input array. Coverages from CoverM contig
    input_options.add_argument(
        '--input',
        help='CoverM coverage results',
        dest="input",
        required=True)

    input_options.add_argument(
        '--assembly',
        help='FASTA file containing scaffolded contigs of the metagenome assembly',
        dest="assembly",
        required=True,
    )

    input_options.add_argument(
        '--min_bin_size',
        help='The minimum size of a returned MAG in base pairs',
        dest="min_bin_size",
        default=200000,
        required=False,
    )

    input_options.add_argument(
        '--scaler',
        help='The method used to scale the input data',
        dest="scaler",
        default='clr',
        choices=['clr', 'minmax', 'none'],
        required=False,
    )

    input_options.add_argument(
        '--min_contig_size',
        help='The minimum contig size to be considered for binning',
        dest="min_contig_size",
        default=1000,
        required=False,
    )

    input_options.add_argument(
        '--output_directory',
        help='Output directory',
        dest="output",
        default="rosella_bins",
        required=False,
    )

    ## UMAP parameters
    input_options.add_argument('--n_neighbors',
                               help='Number of neighbors considered in UMAP',
                               dest="n_neighbors",
                               default=100)

    input_options.add_argument(
        '--min_dist',
        help=
        'Minimum distance used by UMAP during construction of high dimensional graph',
        dest="min_dist",
        default=0)

    input_options.add_argument('--n_components',
                               help='Dimensions to use in UMAP projection',
                               dest="n_components",
                               default=3)

    ## HDBSCAN parameters
    input_options.add_argument('--min_cluster_size',
                               help='Minimum cluster size for HDBSCAN',
                               dest="min_cluster_size",
                               default=5)

    input_options.add_argument('--min_samples',
                               help='Minimum samples for HDBSCAN',
                               dest="min_samples",
                               default=5)

    input_options.add_argument('--cluster_selection_method',
                               help='Cluster selection method used by HDBSCAN. Either "eom" or "leaf"',
                               dest='cluster_selection_method',
                               default='eom')

    ## Genral parameters
    input_options.add_argument(
        '--precomputed',
        help='Minimum cluster size for HDBSCAN',
        dest="precomputed",
        type=str2bool,
        nargs='?',
        const=True,
        default=False,
    )

    input_options.add_argument('--threads',
                               help='Number of threads to run in parallel',
                               dest='threads',
                               default=8)

    ###########################################################################
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parsing input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    if (len(sys.argv) == 2 or len(sys.argv) == 1 or sys.argv[1] == '-h'
            or sys.argv[1] == '--help'):
        phelp()
    else:
        args = main_parser.parse_args()
        time = datetime.datetime.now().strftime('%H:%M:%S %d-%m-%Y')

        if args.log:
            if os.path.isfile(args.log):
                raise Exception("File %s exists" % args.log)
            logging.basicConfig(
                filename=args.log,
                level=debug[args.verbosity],
                format='%(asctime)s %(levelname)s: %(message)s',
                datefmt='%m/%d/%Y %I:%M:%S %p')
        else:
            logging.basicConfig(
                level=debug[args.verbosity],
                format='%(asctime)s %(levelname)s: %(message)s',
                datefmt='%m/%d/%Y %I:%M:%S %p')

        logging.info("Time - %s" % (time))
        logging.info("Command - %s" % ' '.join(sys.argv))

        prefix = args.output
        if not os.path.exists(prefix):
            os.makedirs(prefix)

        if not args.precomputed:
            clusterer = Binner(args.input,
                                prefix,
                                args.assembly,
                                n_neighbors=int(args.n_neighbors),
                                min_cluster_size=int(args.min_cluster_size),
                                min_contig_size=int(args.min_contig_size),
                                min_samples=int(args.min_samples),
                                min_dist=float(args.min_dist),
                                scaler=args.scaler,
                                n_components=int(args.n_components),
                                cluster_selection_method=args.cluster_selection_method,
                                threads=int(args.threads),
                                )
            clusterer.fit_transform()
            clusterer.cluster()
            clusterer.plot()
            clusterer.plot_distances()
            # np.save(prefix + '_labels.npy', clusterer.labels())
            clusterer.bin_contigs(args.assembly, int(args.min_bin_size))
