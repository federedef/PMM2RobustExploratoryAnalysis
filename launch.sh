#!/usr/bin/env bash

CURRENT_PATH=`pwd`
export EXEC_PATH=$FSCRATCH/PMM2REA
export CODE_PATH=$CURRENT_PATH/aux_scripts
export TEMPLATES_PATH=$CURRENT_PATH/templates
export PATH=$CODE_PATH:$PATH
datasets_folder=$CURRENT_PATH/datasets
blacklist="get_synth_data_from_hunter,generate_degs,generate_report_PCA_synth_degs"
blacklist="get_synth_data_from_hunter"

# Parsing inputs

if [ $1 == "parse" ] ; then
    source ~soft_bio_267/initializes/init_python
    source ~/dev_R/ExpHunterSuite/init_degenes_hunter
    if [ ! -d venv ] ; then
            python -m venv venv --system-site-packages
            source ./venv/bin/activate
            pip install -e ~/dev_py/py_exp_calc
            exit 0
    fi
    source ./venv/bin/activate
    
    if [ ! -d datasets ] ; then
        echo -e "Error: You need to create a folder named datasets with"\
                "all datasets needed for the analysis.\n" \
                "These include: phenotypes, miRNA, gene and metabolic data"
    fi
    pushd ./datasets/sample_cohort
    echo -e "Parsing phenotypes\n"
    # head -n 1 phenotypes.txt | tr "\t" "\n"|semtools -i - -O HPO --list_translate names \
    # | cut -f 2 | tr "\n" "\t" | sed -r 's/\t$/\n/g' > phenotypes_n.txt
    # tail -n +2 phenotypes.txt >> phenotypes_n.txt

    echo -e "Parsing genes\n"
    # echo -e "The number of total genes are:"
    # number_of_total_genes=`head -n 1 all_genes_var_filtered.txt | tr -s "\t" "\n" | wc -l`
    # echo -e "$number_of_total_genes"
    # stable_select -i all_genes_var_filtered.txt --label label_genes.txt -o genes.txt
    # number_of_selected_genes=`head -n 1 genes.txt | tr -s "\t" "\n" | wc -l`
    # echo -e "The number of selected genes are: $number_of_selected_genes"
    sed -i "s/^\t//" genes.txt
    head -n 1 genes.txt | tr "\t" "\n" > gene_list.txt
    sed -i '1s/^/ENSEMBL\n/' gene_list.txt
    add_annotation.R -i gene_list.txt -o genes_tr.txt -c 1 -I ENSEMBL -K SYMBOL
    tail -n +2 genes_tr.txt | awk '{if ($2 == "NA") print $1; else print $2 }'  | tr "\n" "\t"| sed -r 's/\t$/\n/g' > genes_n.txt
    tail -n +2 genes.txt >> genes_n.txt

    # head -n 1 miRNA.txt | tr "\t" "\n" > miRNA_list.txt
    # sed -i '1s/^/MIMAT\n/' miRNA_list.txt
    # add_annotation.R -i miRNA_list.txt -o miRNA_tr.txt -c 1 -m
    # tail -n +2 miRNA_tr.txt | awk '{if ($2 == "NA") print $1; else print $2 }'  | tr "\n" "\t"| sed -r 's/\t$/\n/g' > miRNA_n.txt
    # tail -n +2 miRNA.txt >> miRNA_n.txt
    popd

    pushd ./datasets/all_cohort
    head -n 1 phenotypes.txt | tr "\t" "\n"|semtools -i - -O HPO --list_translate names | cut -f 2 | tr "\n" "\t" | sed -r 's/\t$/\n/g' > phenotypes_n.txt
    tail -n +2 phenotypes.txt >> phenotypes_n.txt
    popd
fi

if [ $1 == "ma" ] ; then
	#source ~soft_bio_267/initializes/init_autoflow
	source ~/dev_py/pytoflow/bin/activate
    path_to_autoflow_exec=$EXEC_PATH/multifactor_analysis
    mkdir -p $path_to_autoflow_exec
	variables=`echo -e "
		\\$datasets=$datasets_folder,
        \\$seed=123
	" | tr -d [:space:]`
	AutoFlow -e -w $TEMPLATES_PATH/multifactor_analysis.af -V $variables -o $path_to_autoflow_exec -c 1 -m 60gb -t 0-02:30:00 -n cal $2
fi

if [ $1 == "pa" ] ; then
	#source ~soft_bio_267/initializes/init_autoflow
	source ~/dev_py/pytoflow/bin/activate
    path_to_autoflow_exec=$EXEC_PATH/posterior_analysis
    results_from_ma=$CURRENT_PATH/PCA_top_results #`grep -w "collect_results" $EXEC_PATH/multifactor_analysis/index_execution | cut -f 2`
    mkdir -p $path_to_autoflow_exec
	variables=`echo -e "
		\\$datasets=$datasets_folder,
        \\$results_from_ma=$results_from_ma,
        \\$seed=123
	" | tr -d [:space:]`
	AutoFlow -e -w $TEMPLATES_PATH/posterior_analysis.af -V $variables -o $path_to_autoflow_exec -c 1 -m 20gb -t 0-02:00:00 -n cal $2
fi

if [ "$1" == "check" ] ; then
  	source ~/dev_py/pytoflow/bin/activate
    path_to_autoflow_exec=$EXEC_PATH/multifactor_analysis
	flow_logger -w -e $path_to_autoflow_exec -r all
fi

if [ "$1" == "recover" ] ; then 
    path_to_autoflow_exec=$EXEC_PATH/multifactor_analysis
	flow_logger -w -e $path_to_autoflow_exec --sleep 0.1 -l -p 
fi
