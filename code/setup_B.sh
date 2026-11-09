export R_LIBS_USER=./R_LIBS
mkdir -p R_LIBS
R -e "install.packages(c('igraph', 'optparse'), lib='./R_LIBS', repos='https://cloud.r-project.org')"
