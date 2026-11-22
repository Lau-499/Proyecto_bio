#!/usr/bin/env Rscript
#propiedades de la red y clustering
#lee la red de STRING (02_red_string.R), busca comunidades con Louvain
# y guarda a que cluster pertenece cada gen
#se ejecuta desde la carpeta code: ./03_propiedades_clustering.R

library(optparse)
library(igraph)

option_list <- list(
  make_option(c("-i", "--red"), type = "character", default = "../results/red_string.tsv",
              help = "Archivo con la red (columnas from y to)"),
  make_option(c("-o", "--outdir"), type = "character", default = "../results",
              help = "Carpeta de resultados")
)
opt <- parse_args(OptionParser(option_list = option_list))

#leer la red y pasarla a igraph
red <- read.delim(opt$red, header = TRUE)
g <- graph_from_data_frame(red[, c("from", "to")], directed = FALSE)
cat("Nodos:", gorder(g), "\n")
cat("Aristas:", gsize(g), "\n")

#clustering con Louvain
set.seed(123)
communities <- cluster_louvain(g)
cat("Numero de clusters:", length(communities), "\n")
cat("Modularidad:", modularity(communities), "\n")
print(sizes(communities))

#guardar el cluster de cada gen
clusters <- data.frame(gene = V(g)$name, cluster = membership(communities))
write.table(clusters, file.path(opt$outdir, "clusters.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
