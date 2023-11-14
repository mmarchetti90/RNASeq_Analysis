FROM bioconductor/bioconductor_docker

### INSTALLING PIPELINE PACKAGES ----------- ###

# Installing R packages
RUN Rscript -e 'install.packages(c("BiocManager", \
                                   "ggplot2", \
                                   "msigdbr", \
                                   "pheatmap", \
                                   "RColorBrewer"))'

RUN Rscript -e 'BiocManager::install(c("BiocParallel", \
                                       "clusterProfiler", \
                                       "DESeq2", \
                                       "devtools", \
                                       "DEXSeq", \
                                       "enrichplot", \
                                       "org.Hs.eg.db", \
                                       "org.Mm.eg.db", \
                                       "org.Rn.eg.db", \
                                       "org.Dr.eg.db", \
                                       "org.Dm.eg.db", \
                                       "org.Sc.eg.db", \
                                       "org.Ce.eg.db", \
                                       "org.Ce.eg.db"))'

RUN Rscript -e 'BiocManager::install("pachterlab/sleuth")'

### SETTING WORKING ENVIRONMENT ------------ ###

# Set workdir to /home/
WORKDIR /home/

# Launch bash automatically
CMD ["/bin/bash"]
