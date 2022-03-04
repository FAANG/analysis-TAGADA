FROM nfcore/base:2.1

RUN apt-get update && apt-get install libxt6 ocaml -y

RUN conda update conda -n base -c defaults

RUN conda install mamba -n base -c conda-forge

COPY environment.yml /

RUN mamba env create -f /environment.yml && conda clean -a

RUN rm /environment.yml

ENV PATH /opt/conda/envs/tagada/bin:$PATH

RUN git clone --branch v2.1.7 --depth 1 https://github.com/gpertea/stringtie.git /usr/local/src/StringTie && \
    bash -c 'ln -s /usr/local/src/StringTie/prepDE.py3 /usr/local/bin'

RUN git clone --branch v0.4 --depth 1 https://github.com/sdjebali/Scripts.git /usr/local/src/Scripts && \
    bash -c 'ln -s /usr/local/src/Scripts/* /usr/local/bin'

RUN git clone --branch v1.3 --depth 1 https://github.com/sdjebali/Comptr.git /usr/local/src/Comptr && \
    cd /usr/local/src/Comptr && \
    make && \
    bash -c 'ln -s /usr/local/src/Comptr/comptr /usr/local/bin'

RUN git clone --branch v1.0 --depth 1 https://github.com/sdjebali/Overlap.git /usr/local/src/Overlap && \
    cd /usr/local/src/Overlap && \
    make && \
    bash -c 'ln -s /usr/local/src/Overlap/overlap /usr/local/bin'

RUN git clone -n https://github.com/julienlag/tmerge /usr/local/src/tmerge && \
    cd /usr/local/src/tmerge && \
    git checkout 8b4d6e7c1c94955931946081476e326c4cece161 && \
    bash -c 'ln -s /usr/local/src/tmerge/tmerge /usr/local/bin'

RUN git clone -n https://github.com/julienlag/buildLoci /usr/local/src/buildLoci && \
    cd /usr/local/src/buildLoci && \
    git checkout 0f4ba3c2862b8786a3b94692b877d386d713830b && \
    bash -c 'ln -s /usr/local/src/buildLoci/buildLoci.pl /usr/local/bin'

USER root

RUN git clone --branch v1.0.2 --depth 1 https://github.com/cguyomar/multiqc_feelnc /usr/local/src/multiqc_feelnc && \
    cd /usr/local/src/multiqc_feelnc && \
    python setup.py install

RUN git clone --branch v0.1 --depth 1 https://github.com/cguyomar/custom_images_mqc_plugin /usr/local/src/multiqc_custom_images && \
    cd /usr/local/src/multiqc_custom_images && \
    python setup.py install
