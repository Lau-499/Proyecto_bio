
#!/bin/bash

mkdir -p R_LIBS


R -e "install.packages(c('igraph', 'optparse'), lib='./R_LIBS', repos='https://cloud.r-project.org')"

R -e "install.packages('BiocManager', lib='./R_LIBS', repos='https://cloud.r-project.org')"

R -e ".libPaths(c('./R_LIBS', .libPaths())); BiocManager::install('STRINGdb', lib='./R_LIBS', ask=FALSE, update=FALSE)"
