FROM jupyter/datascience-notebook:r-4.0.3

USER root

COPY package_check.R /home/jovyan/package_check.R

RUN apt-get update && apt-get install -y \
    libglpk-dev \
    libbz2-dev \
    liblzma-dev \
    vim \
 && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('BiocManager', 'kableExtra', 'openxlsx', 'factoextra', 'FactoMineR', 'cowplot', 'ggpubr'), repos = 'https://cran.rstudio.com')"

RUN R -e "BiocManager::install(c('mixOmics', 'PCAtools', 'limma', 'mCSEA'))"
