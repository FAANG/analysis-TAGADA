FROM nfcore/base:1.9

COPY environment.yml /

RUN conda env create -f /environment.yml
RUN conda clean -a
RUN rm /environment.yml

ENV PATH /opt/conda/envs/rnaseq/bin:$PATH

RUN git clone --branch v2.1.4 --depth 1 https://github.com/gpertea/stringtie.git /usr/local/src/StringTie
RUN bash -c 'ln -s /usr/local/src/StringTie/prepDE.py /usr/local/bin'

RUN git clone --depth 1 https://github.com/sdjebali/Scripts.git /usr/local/src/Scripts
RUN bash -c 'ln -s /usr/local/src/Scripts/infer_library_type.sh /usr/local/bin'
RUN bash -c 'ln -s /usr/local/src/Scripts/cutgff.awk /usr/local/bin'
