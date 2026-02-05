source ~soft_bio_267/initializes/init_degenes_hunter

add_annotation.R -i table_genes.txt -o table_genes_n.txt -c 1 -I ENSEMBL -K SYMBOL
awk 'BEGIN{FS=OFS="\t"}{print $5,$2,$3,$4}' table_genes_n.txt > tmp && mv tmp table_genes_n.txt
sed -i '1s/SYMBOL/factor/' table_genes_n.txt
add_annotation.R -i table_mfa.txt -o table_mfa_n.txt -c 1 -I ENSEMBL -K SYMBOL
awk 'BEGIN{FS=OFS="\t"}{print $5,$2,$3,$4}' table_mfa_n.txt > tmp && mv tmp table_mfa_n.txt
sed -i '1s/SYMBOL/factor/' table_mfa_n.txt

