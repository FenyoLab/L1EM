FROM continuumio/miniconda3:4.5.12

WORKDIR /

RUN conda install -y --override-channels -c bioconda -c conda-forge -c defaults python=2.7.15 bwa=0.7.17 samtools=1.9 numpy=1.14.3 scipy=1.1.0 pysam=0.15.0 bedtools=2.27.1
RUN git clone https://github.com/FenyoLab/L1EM/

