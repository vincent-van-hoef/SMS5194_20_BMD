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

# Add Tini. Tini operates as a process subreaper for jupyter. This prevents kernel crashes.
ENV TINI_VERSION v0.6.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]
