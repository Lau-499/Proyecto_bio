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

# 3: Grafica del enriquecimiento global
if (nrow(as.data.frame(enrich_global)) > 0) {
  pdf(file.path(opt$outdir, "GO_global_dotplot.pdf"))
  print(dotplot(enrich_global, showCategory=10, title="GO Enrichment (BP)"))
  dev.off()
}

# 4: Enriquecimiento por cluster
if (!is.null(opt$clusters) && file.exists(opt$clusters)) {
  clusters_table <- read.delim(opt$clusters, header=TRUE, stringsAsFactors=FALSE)
  genes_red <- unique(clusters_table$gene) # Fondo (universe): genes de la red

  genes_clusters <- split(clusters_table$gene, clusters_table$cluster)

  # Filtrar clusters con al menos 5 genes
  min_size <- 5
  genes_clusters <- genes_clusters[lengths(genes_clusters) >= min_size]
  cat("Clusters con al menos", min_size, "genes:", length(genes_clusters), "\n")

  if (length(genes_clusters) > 0) {
    enrichment_results_clusters <- list()

    for (i in names(genes_clusters)) {
      enrichment_results_clusters[[i]] <- enrichGO(
        gene          = genes_clusters[[i]],
        universe      = genes_red,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2
      )

      res_df <- as.data.frame(enrichment_results_clusters[[i]])
      write.csv2(res_df,
                 file=file.path(opt$outdir, paste0("enrichment_cluster_", i, ".csv")),
                 row.names=FALSE)

      if (nrow(res_df) > 0) {
        pdf(file.path(opt$outdir, paste0("dotplot_cluster_", i, ".pdf")), width=9, height=6)
        print(dotplot(enrichment_results_clusters[[i]], showCategory=10,
                      title=paste("GO Enrichment, cluster", i)))
        dev.off()
      }
    }

    # Tabla resumen con el termino mas significativo de cada cluster
    cluster_summary <- data.frame(
      cluster = names(genes_clusters),
      genes   = lengths(genes_clusters),
      top_term = sapply(enrichment_results_clusters, function(x) {
        res <- as.data.frame(x)
        if (nrow(res) == 0) "(none)" else res$Description[1]
      })
    )
    print(cluster_summary)
    write.csv2(cluster_summary,
               file=file.path(opt$outdir, "GO_clusters_summary.csv"),
               row.names=FALSE)
  }
}