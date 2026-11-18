#!/usr/bin/env Rscript
# Script de la red de STRING
# Lee los genes del fenotipo, los busca en STRING y guarda las interacciones.
# Se ejecuta desde la carpeta code: ./02_red_string.R

library(optparse)
library(STRINGdb)

# Opciones del script (como en el ejercicio del manual)
option_list <- list(
  make_option(c("-i", "--genes"), type = "character", default = "../results/genes.tsv",
              help = "Archivo con los genes del HPO"),
  make_option(c("-o", "--outdir"), type = "character", default = "../results",
              help = "Carpeta de resultados"),
  make_option(c("-s", "--score"), type = "integer", default = 800,
              help = "Puntuacion minima de STRING")
)
opt <- parse_args(OptionParser(option_list = option_list))

# 1. Leer los genes
genes <- read.delim(opt$genes, header = TRUE)
genes$name <- trimws(genes$name)
cat("Genes del HPO:", nrow(genes), "\n")

# 2. Buscar los genes en STRING (humano = 9606)
string_db <- STRINGdb$new(version = "12.5", species = 9606,
                          score_threshold = opt$score,
                          input_directory = file.path(opt$outdir, "string_data"))
string_ids <- string_db$map(genes, "name", removeUnmappedRows = TRUE)
cat("\n")

# Genes que no estan en STRING
no_mapeados <- setdiff(genes$name, string_ids$name)
cat("Genes encontrados en STRING:", nrow(string_ids), "\n")
cat("Genes no encontrados:", length(no_mapeados), "\n")
writeLines(no_mapeados, file.path(opt$outdir, "genes_no_mapeados.txt"))

# 3. Interacciones entre nuestros genes
red <- string_db$get_interactions(string_ids$STRING_id)
red <- unique(red[, c("from", "to", "combined_score")])

# Poner el nombre del gen en vez del codigo de STRING
red$from <- string_ids$name[match(red$from, string_ids$STRING_id)]
red$to <- string_ids$name[match(red$to, string_ids$STRING_id)]

# Quitar las repetidas (A-B y B-A es la misma interaccion)
pares <- paste(pmin(red$from, red$to), pmax(red$from, red$to))
red <- red[!duplicated(pares), ]
cat("Interacciones en la red:", nrow(red), "\n")

# Guardar la red
write.table(red, file.path(opt$outdir, "red_string.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
