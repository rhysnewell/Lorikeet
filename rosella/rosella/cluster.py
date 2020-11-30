#!/usr/bin/env python
###############################################################################
# cluster.py - A program which handles the UMAP and HDBSCAN python components
#              of lorikeet
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
import hdbscan
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import skbio.stats.composition
import umap
import numba

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
cluster.py [SUBCOMMAND] ..

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


class Cluster():
    def __init__(
        self,
        count_path,
        output_prefix,
        scaler="clr",
        n_neighbors=20,
        min_dist=0.1,
        n_components=2,
        random_state=42,
        min_cluster_size=100,
        min_samples=50,
        prediction_data=True,
        cluster_selection_method="eom",
        precomputed=False,
        metric='rho',
        hdbscan_metric="euclidean",
    ):

        ## Set up clusterer and UMAP
        self.path = output_prefix
        self.depths = np.load(count_path)

        ## Scale the data
        if scaler.lower() == "minmax":
            self.depths = MinMaxScaler().fit_transform(self.depths)
        elif scaler.lower() == "clr":
            self.depths = skbio.stats.composition.clr(self.depths + 1)
        elif scaler.lower() == "none":
            pass

        if n_neighbors >= int(self.depths.shape[0] * 0.5):
            n_neighbors = max(int(self.depths.shape[0] * 0.5), 2)

        if n_components > self.depths.shape[1]:
            n_components = self.depths.shape[1]

        if metric in ['rho', 'phi', 'phi_dist']:
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

        self.min_samples = min_samples
        self.clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            # min_samples=min_samples,
            prediction_data=prediction_data,
            cluster_selection_method=cluster_selection_method,
            metric=hdbscan_metric,
        )

    def fit_transform(self):
        ## Calculate the UMAP embeddings
        self.embeddings = self.reducer.fit_transform(self.depths)

    def cluster(self):
        ## Cluster on the UMAP embeddings and return soft clusters
        try:
            self.clusterer.fit(self.embeddings)
            self.soft_clusters = hdbscan.all_points_membership_vectors(
                self.clusterer)
            self.soft_clusters = np.array([np.argmax(x) for x in self.soft_clusters])
        except:
            ## Likely integer overflow in HDBSCAN
            ## Try reduce min samples
            self.clusterer = hdbscan.HDBSCAN(
                min_cluster_size=max(int(self.depths.shape[0] * 0.01), 2),
                min_samples=self.min_samples,
                prediction_data=True,
                cluster_selection_method="eom",
            )
            self.clusterer.fit(self.embeddings)
            self.soft_clusters = hdbscan.all_points_membership_vectors(
                self.clusterer)
            self.soft_clusters = np.array([np.argmax(x) for x in self.soft_clusters])


    def cluster_distances(self):
        ## Cluster on the UMAP embeddings and return soft clusters
        try:
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
            self.clusterer.fit(self.depths)

    def plot(self):
        color_palette = sns.color_palette('Paired', 200)
        cluster_colors = [
            color_palette[x] if x >= 0 else (0.5, 0.5, 0.5) for x in self.soft_clusters
        ]
        cluster_member_colors = [
            sns.desaturate(x, p) for x, p in zip(cluster_colors, self.clusterer.probabilities_)
        ]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(self.embeddings[:, 0],
                   self.embeddings[:, 1],
                   s=7,
                   linewidth=0,
                   c=cluster_member_colors,
                   alpha=0.7)
        # ax.add_artist(legend)
        plt.gca().set_aspect('equal', 'datalim')
        plt.title('UMAP projection of variants', fontsize=24)
        plt.savefig(self.path + '_UMAP_projection_with_clusters.png')

    def plot_distances(self):
        self.clusterer.condensed_tree_.plot(
            select_clusters=True,
            selection_palette=sns.color_palette('deep', 20))
        plt.title('Hierarchical tree of clusters', fontsize=24)
        plt.savefig(self.path + '_UMAP_projection_with_clusters.png')

    def labels(self):
        try:
            return self.soft_clusters.astype('int8')
        except AttributeError:
            return self.clusterer.labels_.astype('int8')

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

cluster.py fit --depths depths.npy

''')
    ## Main input array. Depths or Distances
    input_options.add_argument(
        '--input',
        help='.npy file contain depths of variants for each sample',
        dest="input",
        required=True)

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
                               default=1)

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
        prefix = args.input.replace(".npy", "")
        if not args.precomputed:
            clusterer = Cluster(args.input,
                                prefix,
                                n_neighbors=int(args.n_neighbors),
                                min_cluster_size=int(args.min_cluster_size),
                                min_samples=int(args.min_samples),
                                min_dist=float(args.min_dist),
                                n_components=int(args.n_components))
            clusterer.fit_transform()
            clusterer.cluster()
            clusterer.plot()
            np.save(prefix + '_labels.npy', clusterer.labels())
        else:
            clusterer = Cluster(args.input,
                                prefix,
                                n_neighbors=int(args.n_neighbors),
                                min_cluster_size=int(args.min_cluster_size),
                                min_samples=int(args.min_samples),
                                scaler="none",
                                precomputed=args.precomputed)
            clusterer.cluster_distances()
            clusterer.plot_distances()
            np.save(prefix + '_labels.npy', clusterer.labels())
