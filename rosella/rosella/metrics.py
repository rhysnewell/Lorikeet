#!/usr/bin/env python
###############################################################################
# metrics.py - File containing additonal distance metrics for use with UMAP
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
import numba
import numpy as np
import math

###############################################################################                                                                                                                      [44/1010]
################################ - Functions - ################################

@numba.njit()
def tnf(a, b, n_samples):
    # L2 norm is equivalent to euclidean distance
    euc_dist = np.linalg.norm(a[n_samples:] - b[n_samples:])
    return euc_dist

@numba.njit()
def euclidean(a, b, n_samples):
    # Since these compositonal arrays are CLR transformed
    # This is the equivalent to the aitchinson distance but we calculat the l2 norm
    euc_dist = np.linalg.norm(a[:n_samples], b[:n_samples])
    return euc_dist

@numba.njit()
def rho(a, b, n_samples):
    # This is a transformed, inversed version of rho. Normal those -1 <= rho <= 1
    # transformed rho: 0 <= rho <= 2, where 0 is perfect concordance
    covariance_mat = np.cov(a[:n_samples], b[:n_samples], rowvar=True)
    covariance = covariance_mat[0, 1]
    var_a = covariance_mat[0, 0]
    var_b = covariance_mat[1, 1]
    vlr = -2 * covariance + var_a + var_b
    rho = 1 - vlr / (var_a + var_b)
    rho += 1
    rho = 2 - rho
    
    return rho

@numba.njit()
def phi(a, b):
    covariance_mat = np.cov(a, b, rowvar=True)
    covariance = covariance_mat[0, 1]
    var_a = covariance_mat[0, 0]
    var_b = covariance_mat[1, 1]
    phi = 1 + (var_a / var_b) - 2 * np.sqrt(var_a / var_b) * covariance / np.sqrt(var_a * var_b)

    return phi

@numba.njit()
def phi_dist(a, b):
    covariance_mat = np.cov(a, b, rowvar=True)
    covariance = covariance_mat[0, 1]
    var_a = covariance_mat[0, 0]
    var_b = covariance_mat[1, 1]
    phi_dist = abs(math.log(var_a / var_b)) + math.log(2) - math.log(covariance / np.sqrt(var_a * var_b) + 1)

    return phi_dist

@numba.njit()
def aggregate(a, b, n_samples):
    w = n_samples / (n_samples + 1) # weighting by number of samples same as in metabat2
    tnf_dist = tnf(a, b, n_samples)
    aitchinson = euclidean(a, b, n_samples)
    if n_samples >= 3:
        rho_val = rho(a, b, n_samples)
        return tnf_dist^(1-w) * aitchinson^(w) * rho_val
    else:
        return tnf_dist^(1-w) * aitchinson^(w)


