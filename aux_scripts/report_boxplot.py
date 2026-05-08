#!/usr/bin/env python

import argparse
#from py_report_html.py_report_html import Py_report_html 
from py_report_html.py_report_html import Py_report_html
#import Py_report_html
from pathlib import Path

def open_table(file):
    table = []
    with open(file) as f:
        for line in f:
            line = line.strip().split("\t")
            table.append(line)
    return table
# -----------------
parser = argparse.ArgumentParser()
parser.add_argument("--patients2genes", dest= "patients2genes", type=str)
parser.add_argument("--genes", dest= "genes", type=str)
parser.add_argument("--output", dest = "output", type=str)
options = parser.parse_args()
# -----------------

patients2genes = open_table(options.patients2genes)
header_patients2genes = patients2genes.pop(0)
genes = open_table(options.genes)
genes.pop(0)
dim2genes = {}
for row in genes:
    if row[1] in dim2genes:
        dim2genes[row[1]].append(row[0])
    else:
        dim2genes[row[1]] = [row[0]]

final_table = []
for row in patients2genes:
    patient_id = row[0]
    for dim, genes in dim2genes.items():
        for i, gene in enumerate(header_patients2genes):
            if gene in genes:
                final_table.append([patient_id, gene,dim,  patient_id.split("_")[1], row[i+1]])

final_table = sorted(final_table, key= lambda x: x[2])
final_table.insert(0, ["PatID", "Gene", "Dim", "PatGroup", "Counts"])

container = {"final_table": final_table}

report = Py_report_html(container, options.output)
report.compress = True

BASE_DIR = Path(__file__).resolve().parent
file_path = BASE_DIR / "boxplot.txt"
report.build(file_path.read_text())
report.write(options.output + '.html')





