FROM nfcore/base:1.9

RUN apt-get update && apt-get install libxt6 ocaml -y

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN rm /environment.yml

ENV PATH /opt/conda/envs/rnaseq/bin:$PATH

RUN git clone --branch v2.1.4 --depth 1 https://github.com/gpertea/stringtie.git /usr/local/src/StringTie && \
    bash -c 'ln -s /usr/local/src/StringTie/prepDE.py /usr/local/bin'

RUN git clone --branch v0.1 --depth 1 https://github.com/sdjebali/Scripts.git /usr/local/src/Scripts && \
    bash -c 'ln -s /usr/local/src/Scripts/* /usr/local/bin'

RUN git clone --branch v1.3 --depth 1 https://github.com/sdjebali/Comptr.git /usr/local/src/Comptr && \
    cd /usr/local/src/Comptr && \
    make && \
    bash -c 'ln -s /usr/local/src/Comptr/comptr /usr/local/bin'

RUN git clone --branch v1.0 --depth 1 https://github.com/sdjebali/Overlap.git /usr/local/src/Overlap && \
    cd /usr/local/src/Overlap && \
    make && \
    bash -c 'ln -s /usr/local/src/Overlap/overlap /usr/local/bin'

USER root
RUN git clone --branch v1.0.0 --depth 1 https://github.com/cguyomar/multiqc_feelnc /usr/local/src/multiqc_feelnc && \
    cd /usr/local/src/multiqc_feelnc && \
    python setup.py install
