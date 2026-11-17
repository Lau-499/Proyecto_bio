#!/usr/bin/env Rscript
# 02_red_string.R
# Construye la red de interacción de proteínas (STRING) para los genes del
# fenotipo, calcula sus propiedades topológicas y genera las figuras.
#
# Entrada:
#   - Tabla de genes de HPO separada por tabuladores, con columnas "id" y "name"
#     (por defecto ../results/genes.tsv)
#
# Salidas (en la carpeta indicada con -o, por defecto ../results):
#   - genes_no_mapeados.txt     Genes de HPO que no se encuentran en STRING
#
# Estado: versión 1 (lee los genes y los busca en STRING).
# Pendiente: interacciones, propiedades de la red y figuras.
#
# Uso (desde la carpeta code):
#   ./02_red_string.R -i ../results/genes.tsv -o ../results -s 800

suppressPackageStartupMessages({
  library(optparse)
  library(STRINGdb)
})

# Argumentos
option_list <- list(
  make_option(c("-i", "--genes"), type = "character", default = "../results/genes.tsv",
              help = "Tabla de genes de HPO (columnas id y name) [default: %default]"),
  make_option(c("-o", "--outdir"), type = "character", default = "../results",
              help = "Carpeta donde se guardan los resultados [default: %default]"),
  make_option(c("-s", "--score"), type = "integer", default = 800,
              help = "Umbral de puntuación de STRING (0-1000) [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
string_dir <- file.path(opt$outdir, "string_data")
dir.create(string_dir, showWarnings = FALSE)

# 1. Leer los genes del fenotipo
genes <- read.delim(opt$genes, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes$name <- trimws(genes$name)
genes <- genes[genes$name != "", , drop = FALSE]
cat("Genes del fenotipo (HPO):", nrow(genes), "\n")

# 2. Buscar los genes en STRING
# Versión fija para que los resultados no cambien cuando salga una nueva.
string_db <- STRINGdb$new(version = "12.5", species = 9606,
                          score_threshold = opt$score,
                          input_directory = string_dir)

string_ids <- string_db$map(genes, "name", removeUnmappedRows = TRUE)
no_mapeados <- setdiff(genes$name, string_ids$name)
writeLines(no_mapeados, file.path(opt$outdir, "genes_no_mapeados.txt"))
cat("Genes encontrados en STRING:", nrow(string_ids), "\n")
cat("Genes NO encontrados en STRING:", length(no_mapeados), "\n")

