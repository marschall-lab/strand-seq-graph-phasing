Bootstrap: docker
From: continuumio/miniconda3:4.10.3

%files
env_Renv2.yaml /env_Renv2.yaml

%environment
export PATH="/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:"
unset CONDA_DEFAULT_ENV
export ANACONDA_HOME=/opt/conda
export SINGULARITY_SHELL=/bin/bash
export LC_ALL=C

%post
export PATH=/opt/conda/bin:$PATH
conda config --add channels bioconda
conda config --add channels conda-forge
#conda config --add channels r
conda install --yes mamba=0.22.1 pip

mamba env update --quiet -f /env_Renv2.yaml -n base

conda clean --index-cache --tarballs --packages --yes


%labels
AUTHOR Peter Ebert
CONDA v4.10.3
MAMBA v0.22.1
