FROM continuumio/miniconda3:4.12.0

### UPDATING CONDA ------------------------- ###

RUN conda update -y conda

### INSTALLING PIPELINE PACKAGES ----------- ###

# Adding bioconda to the list of channels
RUN conda config --add channels bioconda

# Adding conda-forge to the list of channels
RUN conda config --add channels conda-forge

# Installing mamba
RUN conda install -y mamba

# Installing software
RUN mamba install -y \
    h5py \
    matplotlib \
    numpy \
    pickleshare \
    python-igraph \
    pandas \
    seaborn \
    scipy \
    scikit-learn \
    statsmodels \
    umap-learn && \
    conda clean -afty

RUN mamba install -y -c jmcmurray json

### SETTING WORKING ENVIRONMENT ------------ ###

# Set workdir to /home/
WORKDIR /home/

# Launch bash automatically
CMD ["/bin/bash"]
