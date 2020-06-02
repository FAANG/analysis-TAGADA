FROM nfcore/base:1.9

COPY environment.yml /

RUN conda env create -f /environment.yml
RUN conda clean -a
RUN rm /environment.yml

ENV PATH /opt/conda/envs/rnaseq/bin:$PATH
