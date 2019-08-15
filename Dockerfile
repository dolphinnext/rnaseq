FROM nfcore/base
LABEL author="onur.yukselen@umassmed.edu" description="Docker image containing all requirements for the dolphinnext/rnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
# Install dolphin-tools
RUN mkdir /usr/local/bin/dolphin-tools 
RUN git clone https://github.com/UMMS-biocore/tools /usr/local/bin/dolphin-tools
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext-rnaseq-1.0/bin:/usr/local/bin/dolphin-tools/:$PATH
