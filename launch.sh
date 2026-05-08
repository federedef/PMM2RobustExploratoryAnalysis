#!/usr/bin/env bash

CURRENT_PATH=`pwd`
export EXEC_PATH=$FSCRATCH/PMM2REA
export CODE_PATH=$CURRENT_PATH/aux_scripts
export TEMPLATES_PATH=$CURRENT_PATH/templates
export PATH=$CODE_PATH:$PATH
datasets_folder=$CURRENT_PATH/datasets/datasets
input_folder=$CURRENT_PATH/parsed_datasets
string_cutoff=700
seed=123

if [ $1 == "set_env" ] ; then
    if [ ! -d venv ] ; then
        python -m venv venv --system-site-packages
        source ./venv/bin/activate
        #pip install -e ~/dev_py/py_exp_calc
        pip install 
        echo -e "environment was set up, please run parse option now"
        exit 0
    fi
fi

if [ $1 == "get_string" ] ; then 
    source ~soft_bio_267/initializes/init_python
    wget https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz -O $input_folder/ppi_string.gz
    gzip -d $input_folder/ppi_string.gz
    tr -s " " "\t" < $input_folder/ppi_string | tail -n +2 \
    | awk -v cutoff=$string_cutoff 'BEGIN{FS=OFS="\t"}{if ($3 >= cutoff) print $1,$2}' > tmp && mv tmp $input_folder/ppi_string.txt
    wget https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz -O $input_folder/9606.protein.aliases.v12.0.txt.gz
    gzip -d $input_folder/9606.protein.aliases.v12.0.txt.gz
    grep "Ensembl_gene" $input_folder/9606.protein.aliases.v12.0.txt | cut -f 1,2 > $input_folder/string_ensembl.txt
    grep "Ensembl_HGNC_symbol" $input_folder/9606.protein.aliases.v12.0.txt | cut -f 1,2 > $input_folder/string_hgnc.txt
    rm $input_folder/9606.protein.aliases.v12.0.txt
    cmdtabs -i $input_folder/ppi_string.txt -I $input_folder/string_ensembl.txt -c 1,2 --from 1 --to 2 > $input_folder/ppi_string_ensembl.txt
    cmdtabs -i $input_folder/ppi_string.txt -I $input_folder/string_hgnc.txt -c 1,2 --from 1 --to 2 > $input_folder/ppi_string_hgnc.txt
    rm $input_folder/ppi_string.txt
    rm $input_folder/ppi_string
    rm $input_folder/string_ensembl.txt
    rm $input_folder/string_hgnc.txt
fi

if [ $1 == "parse" ] ; then
    # Parsing inputs
    source ~soft_bio_267/initializes/init_python
    source ~/dev_R/ExpHunterSuite/init_degenes_hunter
    # if [ ! -d venv ] ; then
    #     echo -e "Please run set_env option"
    #     exit 0
    # fi
    # source ./venv/bin/activate
    # if [ ! -d datasets ] ; then
    #     echo -e "Error: You need to create a folder named datasets with"\
    #             "all datasets needed for the analysis.\n" \
    #             "These include: phenotypes, miRNA, gene and metabolic data"
    # fi
    # # Parsing all cohort
    # data_path=$datasets_folder/all_cohort
    # parsed_path=$input_folder/all_cohort
    # mkdir -p $parsed_path

    # ln -s $data_path/phenotypes.txt $parsed_path/phenotypes.txt
    # head -n 1 $parsed_path/phenotypes.txt | tr "\t" "\n"|semtools -i - -O HPO --list_translate names \
    # | cut -f 2 | tr "\n" "\t" | sed -r 's/\t$/\n/g' > $parsed_path/phenotypes_n.txt
    # tail -n +2 $parsed_path/phenotypes.txt >> $parsed_path/phenotypes_n.txt

    # Parsing sample cohort
    data_path=$datasets_folder/sample_cohort
    parsed_path=$input_folder/sample_cohort
    # mkdir -p $parsed_path

    # # Phenotypes
    # echo -e "Parsing phenotypes:\n"
    # ln -s $data_path/phenotypes.txt $parsed_path/phenotypes.txt
    # head -n 1 $parsed_path/phenotypes.txt | tr "\t" "\n"|semtools -i - -O HPO --list_translate names \
    #  | cut -f 2 | tr "\n" "\t" | sed -r 's/\t$/\n/g' > $parsed_path/phenotypes_n.txt
    #  tail -n +2 $parsed_path/phenotypes.txt >> $parsed_path/phenotypes_n.txt

    # Genes
    echo -e "Parsing genes:\n"
    ln -s $data_path/all_genes_var_filtered.txt $parsed_path/all_genes_var_filtered.txt
    head -n 1 $parsed_path/all_genes_var_filtered.txt | tr "\t" "\n" > $parsed_path/all_genes_var_filtered_list.txt
    sed -i '1s/^/ENSEMBL\n/' $parsed_path/all_genes_var_filtered_list.txt
    add_annotation.R -i $parsed_path/all_genes_var_filtered_list.txt -o $parsed_path/all_genes_var_filtered_tr.txt -c 1 -I ENSEMBL -K SYMBOL
    tail -n +2 $parsed_path/all_genes_var_filtered_tr.txt | awk '{if ($2 == "NA") print $1; else print $2 }'  \
    | tr "\n" "\t"| sed -r 's/\t$/\n/g' > $parsed_path/all_genes_var_filtered_n.txt
    tail -n +2 $parsed_path/all_genes_var_filtered.txt >> $parsed_path/all_genes_var_filtered_n.txt
    rm $parsed_path/all_genes_var_filtered_list.txt
    rm $parsed_path/all_genes_var_filtered_tr.txt
    # Gene selection
    # ln -s $data_path/label_genes.txt $parsed_path/label_genes.txt
    # echo -e "The number of total genes are:"
    # number_of_total_genes=`head -n 1 $parsed_path/all_genes_var_filtered.txt | tr -s "\t" "\n" | wc -l`
    # echo -e "$number_of_total_genes"
    # stable_select -i $parsed_path/all_genes_var_filtered.txt --label $parsed_path/label_genes.txt -o $parsed_path/genes.txt
    # number_of_selected_genes=`head -n 1 $parsed_path/genes.txt | tr -s "\t" "\n" | wc -l`
    # echo -e "The number of selected genes are: $number_of_selected_genes"
    # sed -i "s/^\t//" $parsed_path/genes.txt
    # head -n 1 $parsed_path/genes.txt | tr "\t" "\n" > $parsed_path/gene_list.txt
    # sed -i '1s/^/ENSEMBL\n/' $parsed_path/gene_list.txt
    # add_annotation.R -i $parsed_path/gene_list.txt -o $parsed_path/genes_tr.txt -c 1 -I ENSEMBL -K SYMBOL
    # tail -n +2 $parsed_path/genes_tr.txt | awk '{if ($2 == "NA") print $1; else print $2 }'  \
    # | tr "\n" "\t"| sed -r 's/\t$/\n/g' > $parsed_path/genes_n.txt
    # tail -n +2 $parsed_path/genes.txt >> $parsed_path/genes_n.txt

    # miRNA
    echo -e "Parsing miRNAs:\n"
    ln -s $data_path/miRNA.txt $parsed_path/miRNA.txt
    head -n 1 $parsed_path/miRNA.txt | tr "\t" "\n" > $parsed_path/miRNA_list.txt
    sed -i '1s/^/MIMAT\n/' $parsed_path/miRNA_list.txt
    add_annotation.R -i $parsed_path/miRNA_list.txt -o $parsed_path/miRNA_tr.txt -c 1 -m
    tail -n +2 $parsed_path/miRNA_tr.txt | awk '{if ($2 == "NA") print $1; else print $2 }' \
    | tr "\n" "\t"| sed -r 's/\t$/\n/g' > $parsed_path/miRNA_n.txt
    tail -n +2 $parsed_path/miRNA.txt >> $parsed_path/miRNA_n.txt

    # # metabolomics
    # ln -s $data_path/metabolomics.txt $parsed_path/metabolomics.txt

    # # Supp
    # for file in severity_scales.txt severity.txt variants.txt ; do 
    #     ln -s $data_path/$file $parsed_path/$file
    # done
fi

if [ $1 == "ma" ] ; then
	#source ~soft_bio_267/initializes/init_autoflow
	source ~soft_bio_267/initializes/init_python
    path_to_autoflow_exec=$EXEC_PATH/multifactor_analysis
    mkdir -p $path_to_autoflow_exec
	variables=`echo -e "
		\\$datasets=$input_folder,
        \\$seed=$seed
	" | tr -d [:space:]`
	AutoFlow -e -w $TEMPLATES_PATH/multifactor_analysis.af -V $variables -o $path_to_autoflow_exec -c 1 -m 60gb -t 0-02:30:00 -n cal $2
fi

if [ $1 == "pa" ] ; then
	#source ~soft_bio_267/initializes/init_autoflow
	source ~soft_bio_267/initializes/init_python
    path_to_autoflow_exec=$EXEC_PATH/posterior_analysis
    results_from_ma=$CURRENT_PATH/PCA_top_results #`grep -w "collect_results" $EXEC_PATH/multifactor_analysis/index_execution | cut -f 2`
    mkdir -p $path_to_autoflow_exec
	variables=`echo -e "
		\\$datasets=$input_folder,
        \\$results_from_ma=$results_from_ma,
        \\$correlation_cutoff=0.70,
        \\$seed=$seed,
        \\$database_mirna=mirecords
	" | tr -d [:space:]`
	AutoFlow -e -w $TEMPLATES_PATH/posterior_analysis.af -V $variables -o $path_to_autoflow_exec -c 1 -m 20gb -t 0-02:00:00 -n cal $2
fi

if [ "$1" == "check" ] ; then
  	source ~soft_bio_267/initializes/init_python
    path_to_autoflow_exec=$EXEC_PATH/multifactor_analysis
	flow_logger -w -e $path_to_autoflow_exec -r all
fi

if [ "$1" == "recover" ] ; then
    source ~soft_bio_267/initializes/init_python
    path_to_autoflow_exec=$EXEC_PATH/multifactor_analysis
	flow_logger -w -e $path_to_autoflow_exec --sleep 0.1 -l -p 
fi
