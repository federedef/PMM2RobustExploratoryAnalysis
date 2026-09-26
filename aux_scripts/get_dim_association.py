#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd

from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr


def load_pc_table(path):
    """
    Load PC/loadings table.

    Expected format:

                Dim.1    Dim.2    Dim.3 ...
    feature1    ...
    feature2    ...

    The first column is used as the feature identifier.
    """

    df = pd.read_csv(
        path,
        sep="\t",
        index_col=0
    )

    pc_columns = [
        col for col in df.columns
        if str(col).startswith("Dim.")
    ]

    if not pc_columns:
        raise ValueError(
            f"No Dim.* columns detected in {path}"
        )

    return df[pc_columns].astype(float)


def get_pc_correlations(reference, comparison):
    """
    Compute Pearson correlation between every reference PC
    and every comparison PC using common features.
    """

    common_features = reference.index.intersection(
        comparison.index
    )

    if len(common_features) < 3:
        raise ValueError(
            "Fewer than 3 common features between tables."
        )

    ref = reference.loc[common_features]
    comp = comparison.loc[common_features]

    corr = pd.DataFrame(
        index=ref.columns,
        columns=comp.columns,
        dtype=float
    )

    for ref_pc in ref.columns:
        for comp_pc in comp.columns:

            x = ref[ref_pc].values
            y = comp[comp_pc].values

            r, _ = pearsonr(x, y)

            corr.loc[ref_pc, comp_pc] = r

    return corr, len(common_features)


def match_pcs(corr):
    """
    Find one-to-one PC correspondence maximizing absolute
    Pearson correlation.

    Example:
        Reference Dim.1 -> LOO Dim.2
        Reference Dim.2 -> LOO Dim.1
    """

    abs_corr = np.abs(corr.values)

    # Hungarian algorithm minimizes cost,
    # so convert maximum correlation into minimum cost.
    cost = -abs_corr

    ref_idx, comp_idx = linear_sum_assignment(cost)

    matches = {}

    for i, j in zip(ref_idx, comp_idx):

        ref_pc = corr.index[i]
        comp_pc = corr.columns[j]

        matches[ref_pc] = {
            "matched_pc": comp_pc,
            "correlation": corr.iloc[i, j],
            "abs_correlation": abs(corr.iloc[i, j])
        }

    return matches


parser = argparse.ArgumentParser(
    description=(
        "Match PCs from multiple analyses to PCs from "
        "a reference analysis using loadings correlations."
    )
)

parser.add_argument(
    "--input_reference",
    type=str,
    required=True,
    help="Reference PC/loadings table."
)

parser.add_argument(
    "--input_list",
    type=str,
    required=True,
    help=(
        "Comma-separated list of PC/loadings tables "
        "to compare against the reference."
    )
)

parser.add_argument(
    "--output",
    type=str,
    required=True,
    help="Output table containing PC correspondences."
)

parser.add_argument(
    "--correlation_output",
    type=str,
    default=None,
    help=(
        "Optional output table containing the correlation "
        "of every matched PC."
    )
)

options = parser.parse_args()


# -------------------------
# Load reference
# -------------------------

reference = load_pc_table(
    options.input_reference
)

reference_pcs = list(reference.columns)


# -------------------------
# Output structures
# -------------------------

mapping_rows = []
correlation_rows = []


# Reference row
reference_row = {
    "file": "reference"
}

for pc in reference_pcs:
    reference_row[pc] = pc

mapping_rows.append(reference_row)


# -------------------------
# Compare each analysis
# -------------------------

for input_file in options.input_list.split(","):

    comparison = load_pc_table(
        input_file
    )

    corr_matrix, n_common = get_pc_correlations(
        reference,
        comparison
    )

    matches = match_pcs(
        corr_matrix
    )

    name = os.path.basename(input_file)

    mapping_row = {
        "file": name
    }

    correlation_row = {
        "file": name,
        "common_features": n_common
    }

    for ref_pc in reference_pcs:

        if ref_pc in matches:

            match = matches[ref_pc]

            mapping_row[ref_pc] = match["matched_pc"]

            correlation_row[
                ref_pc
            ] = match["correlation"]

        else:

            mapping_row[ref_pc] = "NA"
            correlation_row[ref_pc] = np.nan

    mapping_rows.append(mapping_row)
    correlation_rows.append(correlation_row)


# -------------------------
# Write PC mapping
# -------------------------

mapping_df = pd.DataFrame(
    mapping_rows
)

mapping_df.to_csv(
    options.output,
    sep="\t",
    index=False
)


# -------------------------
# Optional correlation file
# -------------------------

if options.correlation_output is not None:

    correlation_df = pd.DataFrame(
        correlation_rows
    )

    correlation_df.to_csv(
        options.correlation_output,
        sep="\t",
        index=False
    )


# -------------------------
# Print
# -------------------------

print("\nPC correspondence:")
print(mapping_df.to_string(index=False))

if options.correlation_output is not None:

    print("\nMatched-PC correlations:")
    print(correlation_df.to_string(index=False))