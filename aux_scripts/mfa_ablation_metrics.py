#!/usr/bin/env python

import numpy as np
import argparse
from scipy.spatial.distance import pdist, squareform

def rv_coefficient(X, Y):
    """
    Classical Escoufier RV coefficient.

    X and Y:
        rows = same samples
        columns = dimensions/features

    X and Y may contain different numbers of columns.
    """

    # Column-center both configurations
    X = X - X.mean(axis=0, keepdims=True)
    Y = Y - Y.mean(axis=0, keepdims=True)

    XX = X @ X.T
    YY = Y @ Y.T

    numerator = np.trace(XX @ YY)
    denominator = np.sqrt(
        np.trace(XX @ XX) *
        np.trace(YY @ YY)
    )

    return numerator / denominator


def rv2_coefficient(X, Y):
    """
    Modified RV coefficient (RV2).

    Removes diagonal elements from the sample cross-product matrices.
    This can be useful with very small n because classical RV can be
    upward biased.
    """

    X = X - X.mean(axis=0, keepdims=True)
    Y = Y - Y.mean(axis=0, keepdims=True)

    XX = X @ X.T
    YY = Y @ Y.T

    np.fill_diagonal(XX, 0.0)
    np.fill_diagonal(YY, 0.0)

    numerator = np.sum(XX * YY)

    denominator = np.sqrt(
        np.sum(XX ** 2) *
        np.sum(YY ** 2)
    )

    return numerator / denominator

def permanova_r2(coords, groups):
    """
    Calculate PERMANOVA R2 only.

    Parameters
    ----------
    coords : array-like, shape (n_samples, n_dimensions)
        Sample coordinates in MFA space.

    groups : array-like, shape (n_samples,)
        Group labels, e.g. H / L.

    Returns
    -------
    r2 : float
        Proportion of multivariate variation associated with group membership.
    """

    coords = np.asarray(coords, dtype=float)
    groups = np.asarray(groups)

    # Euclidean distance matrix between samples
    D = squareform(pdist(coords, metric="euclidean"))

    n = len(groups)

    # Total sum of squares
    upper = np.triu_indices(n, k=1)
    ss_total = np.sum(D[upper] ** 2) / n

    # Within-group sum of squares
    ss_within = 0.0

    for group in np.unique(groups):

        idx = np.where(groups == group)[0]
        ng = len(idx)

        sub_D = D[np.ix_(idx, idx)]
        upper_group = np.triu_indices(ng, k=1)

        ss_within += (
            np.sum(sub_D[upper_group] ** 2) / ng
        )

    # Between-group variability
    ss_between = ss_total - ss_within

    # PERMANOVA effect size
    r2 = ss_between / ss_total

    return r2

def load_coords(file):
    smp2coords = {}
    with open(file, "r") as f:
        n_dims = 2
        for index, line in enumerate(f):
            line = line.strip().split("\t")
            if index == 0:
                n_dims = len(line) - 2
                continue
            smp_id = line[-1]
            dims = [float(coord) for coord in line[0:n_dims]]
            smp2coords[smp_id] = dims
    return smp2coords

def load_groups(file):
    smp2groups = {}
    with open(file, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            smp2groups[line[0]] = line[1]
    return smp2groups

# options
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_full", dest="input_full", type=str, required=True)
parser.add_argument("-j", "--input_ablation", dest="input_ablation", type=str, required=True)
parser.add_argument("-g", "--groups", dest="groups", type=str, required=True)
options = parser.parse_args()

# main
coords_ablation = load_coords(options.input_ablation)
coords_full = load_coords(options.input_full)

# Creating corresponding matrices:
ids = list(coords_ablation.keys())
X = np.array([coords_ablation[smp] for smp in ids])
Y = np.array([coords_full[smp] for smp in ids])
groups = load_groups(options.groups)
G = np.array([groups[smp] for smp in ids])

# metrics
rv = rv_coefficient(X, Y)
rv2 = rv2_coefficient(X, Y)
print(f"RV coefficient:\t{np.round(rv, 3)}")
print(f"RV2 coefficient:\t{np.round(rv2, 3)}")
permanova_r2_metric = permanova_r2(X, G)
permanova_r2_full = permanova_r2(Y, G)
r2_diff = permanova_r2_metric - permanova_r2_full
print(f"R2 difference:\t{np.abs(np.round(r2_diff, 3))}")

# print("the metrics are:", rv, rv2, permanova_r2)



