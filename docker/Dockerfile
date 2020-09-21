FROM continuumio/miniconda3

MAINTAINER Lijun Yao < lijunyao@wustl.edu>

COPY scrna_env.yml ./
COPY myscript.sh ./
COPY perl ./

RUN conda create -n scrna python=3.6
RUN echo "source activate scrna" > ~/.bashrc
ENV PATH /opt/conda/envs/scrna/bin:$PATH

#RUN conda env create -f scrna_env.yml
##RUN echo "source activate scrna" > ~/.bashrc
#ENV PATH /opt/conda/envs/scrna/bin:$PATH

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install -c bioconda r-seurat
RUN conda install -c bioconda r-optparse

RUN conda install -c r r-getopt
RUN conda install -c r r-dplyr
RUN conda install -c bioconda bioconductor-biomaRt

#RUN conda install pandas
#RUN conda install numpy
#RUN conda install scipy
RUN pip install numpy
RUN pip install pandas 
RUN pip install scipy
RUN pip install statsmodels
RUN conda install -c bioconda perl-statistics-distributions
RUN conda install -c bioconda perl-statistics-descriptive
RUN conda install perl-app-cpanminus
RUN cpanm Math::SigFigs

WORKDIR /opt
COPY downloaded_db ./DOWNLOAD_DB
COPY GTEX_lib ./GTEX_lib
#COPY myscript.sh ./
COPY tumor_specific_markers_V1.0.R ./
COPY step0.get_surface_genes_metadata.pl ./
COPY step1.get_surface_genes_expression_data.pl ./
COPY step2.analyze_expression_change.pl ./
COPY HPA_Protein_filtering.py ./
COPY HPA_RNA_find_outlier_expressions_t_test.pl ./
COPY adding_tissue_specificity_annotations.R ./


#CMD /bin/bash
