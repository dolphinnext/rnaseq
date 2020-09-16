# Don't upgrade nfcore/base, it creates "Kernel too old" error for singularity (because of the debian image)
FROM nfcore/base:1.7 

LABEL author="onur.yukselen@umassmed.edu" description="Docker image containing all requirements for the dolphinnext/rnaseq pipeline"

RUN apt-get update && apt-get install -y gcc 
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
# Install dolphin-tools
RUN git clone https://github.com/dolphinnext/tools /usr/local/bin/dolphin-tools
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext-rnaseq-2.0/bin:/usr/local/bin/dolphin-tools/:$PATH

# Install tophat-2.1.1
RUN wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz && tar -xvzf tophat-2.1.1.Linux_x86_64.tar.gz && mv tophat-2.1.1.Linux_x86_64/ /usr/local/bin/dolphin-tools/tophat-2.1.1/
export PATH=/usr/local/bin/dolphin-tools/tophat-2.1.1:$PATH