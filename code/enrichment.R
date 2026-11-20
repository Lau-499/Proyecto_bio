library(optparse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

option_list <- list(
  make_option(c("-g", "--genes"), type="character", default=NULL,
              help="Fichero con la lista de genes (genes.tsv)"),
  make_option(c("-c", "--clusters"), type="character", default=NULL,
              help="Fichero clusters.tsv (gen y cluster). Opcional"),
  make_option(c("-o", "--outdir"), type="character", default="../results/enrichment",
              help="Carpeta de salida")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# 1: Lectura de genes
genes_table <- read.table(opt$genes, sep="\t", header=TRUE, strip.white=TRUE)
genes <- trimws(genes_table[, ncol(genes_table)])
cat("Numero de genes:", length(genes), "\n")

# 2: Enriquecimiento global
enrich_global <- enrichGO(gene=genes, OrgDb=org.Hs.eg.db, keyType="SYMBOL",
                          ont="BP", pAdjustMethod="BH",
                          pvalueCutoff=0.05, qvalueCutoff=0.2)
write.csv2(as.data.frame(enrich_global),
           file=file.path(opt$outdir, "GO_global.csv"), row.names=FALSE)