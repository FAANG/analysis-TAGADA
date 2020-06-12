FROM nfcore/base:1.9

COPY environment.yml /

RUN conda env create -f /environment.yml
RUN conda clean -a
RUN rm /environment.yml

ENV PATH /opt/conda/envs/rnaseq/bin:$PATH

RUN git clone --depth=1 https://github.com/sdjebali/Scripts.git /usr/local/src/scripts
RUN bash -c 'ln -s /usr/local/src/scripts/*/*.{awk,py,sh,R} /usr/local/bin'
