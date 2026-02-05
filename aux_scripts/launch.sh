#!/usr/bin/env bash
#source ~soft_bio_267/initializes/init_python
#source ./venv/bin/activate
source ./venv_new/bin/activate
./report_boxplot.py --patients2genes "mirna_n.txt" --genes "table_mirna_n.txt" --output "mirna"
./report_boxplot.py --patients2genes "genes_n.txt" --genes "table_genes_n.txt" --output "genes"
./report_boxplot.py --patients2genes "genes_mirna_n.txt" --genes "table_mfa_n.txt" --output "mfa"

