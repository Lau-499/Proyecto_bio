#!/bin/bash
cd "$(dirname "$0")"

mkdir -p "$PWD/R_LIBS"
echo "CXX11STD = -std=gnu++17" > "$PWD/Makevars_local"
export R_LIBS_USER="$PWD/R_LIBS"
export R_MAKEVARS_USER="$PWD/Makevars_local"
export USE_BUNDLED_LIBUV=1

# Paquetes base y red
R -e "install.packages(c('igraph', 'optparse', 'BiocManager'), lib='$PWD/R_LIBS', repos='https://cloud.r-project.org')"
R -e "BiocManager::install('STRINGdb', lib='$PWD/R_LIBS', ask=FALSE, update=FALSE)"

# Paquetes de enriquecimiento (enrichment.R)
R -e "options(repos=c(CRAN='https://packagemanager.posit.co/cran/2024-03-01')); install.packages('https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.5.2.tar.gz', repos=NULL, type='source', lib='$PWD/R_LIBS')"
R -e "options(repos=c(CRAN='https://packagemanager.posit.co/cran/2024-03-01')); BiocManager::install(c('org.Hs.eg.db', 'ggtree', 'enrichplot', 'clusterProfiler'), lib='$PWD/R_LIBS', update=FALSE, ask=FALSE)"

rm -f "$PWD/Makevars_local"