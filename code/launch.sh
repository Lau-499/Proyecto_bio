
#!/bin/bash

export R_LIBS_USER=./R_LIBS

mkdir -p ../results

mkdir -p ../results/string_data

# Descargar los genes asociados al embarazo ectopico
wget "https://ontology.jax.org/api/network/annotation/HP:0031456/download/gene" -O ../results/hpo_ectopic.txt

# Quitar espacios del archivo
tr -d ' ' < ../results/hpo_ectopic.txt > ../results/genes.tsv

# Ejecutar el analisis de STRING
Rscript 02_red_string.R \
-i ../results/genes.tsv \
-o ../results \
-s 800

# Propiedades de la red y clustering
Rscript 03_propiedades_clustering.R \
-i ../results/red_string.tsv \
-o ../results

# Ejecutar el analisis de enriquecimiento (global y por clusters)
Rscript enrichment.R \
-g ../results/genes.tsv \
-c ../results/clusters.tsv \
-o ../results/enrichment