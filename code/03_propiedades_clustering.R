#!/usr/bin/env Rscript
# Propiedades de la red y clustering
# Lee la red de STRING (02_red_string.R), calcula sus propiedades,
# busca comunidades con Louvain y guarda a que cluster pertenece cada gen.
# Se ejecuta desde la carpeta code: ./03_propiedades_clustering.R

library(optparse)
library(igraph)

option_list <- list(
  make_option(c("-i", "--red"), type = "character", default = "../results/red_string.tsv",
              help = "Archivo con la red (columnas from y to)"),
  make_option(c("-o", "--outdir"), type = "character", default = "../results",
              help = "Carpeta de resultados")
)
opt <- parse_args(OptionParser(option_list = option_list))

# 1. Leer la red y pasarla a igraph
red <- read.delim(opt$red, header = TRUE)
g <- graph_from_data_frame(red[, c("from", "to")], directed = FALSE)
cat("Nodos:", gorder(g), "\n")
cat("Aristas:", gsize(g), "\n")

# 2. Propiedades de la red
grados <- degree(g)
densidad <- edge_density(g)
propiedades <- data.frame(
  propiedad = c("nodos", "aristas", "grado_medio", "densidad", "dispersion",
                "componentes", "transitividad", "longitud_media_camino",
                "diametro", "asortatividad"),
  valor = c(gorder(g), gsize(g), mean(grados), densidad, 1 - densidad,
            components(g)$no, transitivity(g), mean_distance(g),
            diameter(g), assortativity_degree(g))
)
print(propiedades)
write.table(propiedades, file.path(opt$outdir, "propiedades_red.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Centralidad de cada gen
# La red no es conexa, asi que usamos la centralidad armonica en vez de closeness
centralidad <- data.frame(gene = V(g)$name, grado = grados,
                          armonica = harmonic_centrality(g))
centralidad <- centralidad[order(centralidad$grado, decreasing = TRUE), ]
cat("Genes con mas conexiones:\n")
print(head(centralidad, 10))
write.table(centralidad, file.path(opt$outdir, "centralidad_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 3. Clustering con Louvain
set.seed(123)
communities <- cluster_louvain(g)
cat("Numero de clusters:", length(communities), "\n")
cat("Modularidad:", modularity(communities), "\n")
print(sizes(communities))

# 4. Guardar el cluster de cada gen
clusters <- data.frame(gene = V(g)$name, cluster = membership(communities))
write.table(clusters, file.path(opt$outdir, "clusters.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 5. Figura de los clusters
set.seed(123)
layout <- layout_with_fr(g)

pdf(file.path(opt$outdir, "red_clusters.pdf"))
plot(communities, g, layout = layout, vertex.size = 8,
     vertex.label.cex = 0.6, main = "Clustering con Louvain")
dev.off()

# 6. Comparar la modularidad con redes aleatorias
set.seed(123)
q_random <- replicate(100, modularity(cluster_louvain(sample_gnm(gorder(g), gsize(g)))))
q_red <- modularity(communities)
cat("Modularidad media en redes aleatorias:", mean(q_random), "\n")
cat("Redes aleatorias con Q mayor o igual:", sum(q_random >= q_red), "de 100\n")

pdf(file.path(opt$outdir, "modularidad_aleatoria.pdf"))
hist(q_random, xlim = c(0, 0.5), col = "grey80", main = "Q en redes aleatorias equivalentes",
     xlab = "Modularidad Q")
abline(v = q_red, col = "red", lwd = 3)
dev.off()
