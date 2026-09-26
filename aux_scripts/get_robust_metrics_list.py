#!/usr/bin/env python3

import argparse
from py_cmdtabs.cmdtabs import CmdTabs
import numpy as np

def obtain_jaccard(set_a, set_b):
    union = set_a | set_b
    if len(union) == 0:
        return 1.0
    return len(set_a & set_b) / len(union)

def obtain_retention(reference_set, test_set):
    return len(reference_set & test_set) / len(reference_set)


def obtain_simpson(set_a, set_b):
    min_size = min(len(set_a), len(set_b))

    if min_size == 0:
        if len(set_a) == 0 and len(set_b) == 0:
            return 1.0
        return 0.0

    return len(set_a & set_b) / min_size

# options
parser = argparse.ArgumentParser(description='Get the robust metrics list.')
parser.add_argument('--input_list', type=str, required=True, help='Comma-separated list of input files containing the top genes from different LOO analyses.')
parser.add_argument('--input_reference', type=str, required=True, help='Input file containing the reference list of top genes.')
options = parser.parse_args()

# main
ref_set = set(el[0] for el in CmdTabs.load_input_data(options.input_reference))
metrics = {"jaccards":[], "retentions":[], "simpsons":[], "num_of_elements": []}
for input_file in options.input_list.split(','):
    elements = CmdTabs.load_input_data(input_file)
    elements = set(el[0] for el in elements)
    metrics["jaccards"].append(obtain_jaccard(ref_set, elements))
    metrics["retentions"].append(obtain_retention(ref_set, elements))
    metrics["simpsons"].append(obtain_simpson(ref_set, elements))
    metrics["num_of_elements"].append(len(elements))


print("number of elements in ref set", len(ref_set))
print("number of elements in elements", np.median(metrics["num_of_elements"]))

print("jaccard mean", np.mean(metrics["jaccards"]))
print("retention mean", np.mean(metrics["retentions"]))
print("simpson mean", np.mean(metrics["simpsons"]))

print("jaccard median", np.median(metrics["jaccards"]))
print("retention median", np.median(metrics["retentions"]))
print("simpson median", np.median(metrics["simpsons"]))




