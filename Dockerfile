FROM nfcore/base:1.9

RUN apt-get update && apt-get install libxt6

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN rm /environment.yml

ENV PATH /opt/conda/envs/rnaseq/bin:$PATH

RUN git clone --branch v2.1.4 --depth 1 https://github.com/gpertea/stringtie.git /usr/local/src/StringTie && \
    bash -c 'ln -s /usr/local/src/StringTie/prepDE.py /usr/local/bin'

RUN git clone --depth 1 https://github.com/sdjebali/Scripts.git /usr/local/src/Scripts && \
    bash -c 'ln -s /usr/local/src/Scripts/infer_library_type.sh /usr/local/bin' && \
    bash -c 'ln -s /usr/local/src/Scripts/cutgff.awk /usr/local/bin'
