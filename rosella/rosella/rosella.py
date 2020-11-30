#!/usr/bin/env python
###############################################################################
# rosella.py - A fast binning algorithm spinning off of the methodology of
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
import datetime

# Function imports
import numpy as np

# Self imports
from .binning import Binner
from .cluster import Cluster

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
def main():
    ############################ ~ Main Parser ~ ##############################
    main_parser = argparse.ArgumentParser(prog='rosella',
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
    fit_options = subparsers.add_parser(
        'fit',
        description='Perform UMAP and then HDBSCAN on array of variant depths',
        formatter_class=CustomHelpFormatter,
        epilog='''
                                    ~ fit ~
        How to use fit:

        rosella fit --depths depths.npy

        ''')
    ## Main input array. Depths or Distances
    fit_options.add_argument(
        '--input',
        help='.npy file contain depths of variants for each sample',
        dest="input",
        required=True)

    ## UMAP parameters
    fit_options.add_argument('--n_neighbors',
                             help='Number of neighbors considered in UMAP',
                             dest="n_neighbors",
                             default=100)

    fit_options.add_argument(
        '--min_dist',
        help=
        'Minimum distance used by UMAP during construction of high dimensional graph',
        dest="min_dist",
        default=0)

    fit_options.add_argument('--n_components',
                             help='Dimensions to use in UMAP projection',
                             dest="n_components",
                             default=3)

    fit_options.add_argument('--metric',
                             help='Metric to use in UMAP projection',
                             dest="metric",
                             default="rho")

    ## HDBSCAN parameters
    fit_options.add_argument('--min_cluster_size',
                             help='Minimum cluster size for HDBSCAN',
                             dest="min_cluster_size",
                             default=5)

    fit_options.add_argument('--min_samples',
                             help='Minimum samples for HDBSCAN',
                             dest="min_samples",
                             default=1)

    ## Genral parameters
    fit_options.add_argument(
        '--precomputed',
        help='Minimum cluster size for HDBSCAN',
        dest="precomputed",
        type=str2bool,
        nargs='?',
        const=True,
        default=False,
    )

    fit_options.set_defaults(func=fit)


    bin_options = subparsers.add_parser(
        'bin',
        description='Perform UMAP and then HDBSCAN on array of variant depths',
        formatter_class=CustomHelpFormatter,
        epilog='''
                                ~ bin ~
    How to use bin:

    rosella bin --input coverm_output.tsv --assembly scaffolds.fasta

    ''')
    ## Main input array. Coverages from CoverM contig
    bin_options.add_argument(
        '--input',
        help='CoverM coverage results',
        dest="input",
        required=True)

    bin_options.add_argument(
        '--assembly',
        help='FASTA file containing scaffolded contigs of the metagenome assembly',
        dest="assembly",
        required=True,
    )

    bin_options.add_argument(
        '--min_bin_size',
        help='The minimum size of a returned MAG in base pairs',
        dest="min_bin_size",
        default=200000,
        required=False,
    )

    bin_options.add_argument(
        '--scaler',
        help='The method used to scale the input data',
        dest="scaler",
        default='clr',
        choices=['clr', 'minmax', 'none'],
        required=False,
    )

    bin_options.add_argument(
        '--min_contig_size',
        help='The minimum contig size to be considered for binning',
        dest="min_contig_size",
        default=1000,
        required=False,
    )

    bin_options.add_argument(
        '--output_directory',
        help='Output directory',
        dest="output",
        default="rosella_bins",
        required=False,
    )

    ## UMAP parameters
    bin_options.add_argument('--n_neighbors',
                             help='Number of neighbors considered in UMAP',
                             dest="n_neighbors",
                             default=100)

    bin_options.add_argument(
        '--min_dist',
        help=
        'Minimum distance used by UMAP during construction of high dimensional graph',
        dest="min_dist",
        default=0)

    bin_options.add_argument('--n_components',
                             help='Dimensions to use in UMAP projection',
                             dest="n_components",
                             default=3)

    bin_options.add_argument('--metric',
                             help='Metric to use in UMAP projection',
                             dest="metric",
                             default="aggregate")
    ## HDBSCAN parameters
    bin_options.add_argument('--min_cluster_size',
                             help='Minimum cluster size for HDBSCAN',
                             dest="min_cluster_size",
                             default=5)

    bin_options.add_argument('--min_samples',
                             help='Minimum samples for HDBSCAN',
                             dest="min_samples",
                             default=5)

    bin_options.add_argument('--cluster_selection_method',
                             help='Cluster selection method used by HDBSCAN. Either "eom" or "leaf"',
                             dest='cluster_selection_method',
                             default='eom')

    ## Genral parameters
    bin_options.add_argument(
        '--precomputed',
        help='Minimum cluster size for HDBSCAN',
        dest="precomputed",
        type=str2bool,
        nargs='?',
        const=True,
        default=False,
    )

    bin_options.add_argument('--threads',
                             help='Number of threads to run in parallel',
                             dest='threads',
                             default=8)
    bin_options.set_defaults(func=bin)

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

        args.func(args)

def fit(args):
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


def bin(args):
    prefix = args.output
    if not os.path.exists(prefix):
        os.makedirs(prefix)

    if not args.precomputed:
        clusterer = Binner(args.input,
                           prefix,
                           args.assembly,
                           n_neighbors=int(args.n_neighbors),
                           metric=args.metric,
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

def phelp():
    print("""
Usage:
rosella [SUBCOMMAND] ..

Subcommands:
bin - Bin sets of metagenomic contigs into MAGs
fit - Genotype variants into metagenomic strains *For use with Lorikeet*

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



if __name__ == '__main__':

   sys.exit(main())
