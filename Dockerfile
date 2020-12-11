FROM rocker/verse:4.0.3

USER root

COPY package_check.R /package_check.R

RUN apt-get update && apt-get install -y \
    libglpk-dev \
    libbz2-dev \
    liblzma-dev \
    python3.8 \
    python3-pip \
    python3-setuptools \
    python3-dev \
 && rm -rf /var/lib/apt/lists/*

RUN pip3 install -U jupyter-book ghp-import jupytext

RUN R -e "install.packages(c('devtools', 'kableExtra', 'openxlsx', 'factoextra', 'FactoMineR', 'cowplot', 'ggpubr'), repos = 'https://cran.rstudio.com')"

RUN R -e "BiocManager::install(c('mixOmics', 'PCAtools', 'limma', 'mCSEA'))"

RUN R -e "devtools::install_github('IRkernel/IRkernel'); IRkernel::installspec()"