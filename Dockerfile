FROM ensemblorg/ensembl-vep

USER root

ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

# Updates
RUN apt-get update
RUN apt-get install -y --no-install-recommends apt-utils
RUN apt-get -y upgrade
RUN apt-get install -y wget build-essential libz-dev libcurl3-dev gcc libssl-dev libbz2-dev

# conda installation
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  -O ./miniconda.sh && \
    /bin/bash ./miniconda.sh -b -p /opt/conda && \
    rm ./miniconda.sh

# python pip updates
RUN conda install python=3.6
RUN conda config --add channels bioconda
RUN conda install cyvcf2 -y
RUN conda install cython -y
RUN conda install pybigwig -y
RUN pip install --upgrade pip
RUN python -V
RUN pip -V

COPY . ./mmsplice
WORKDIR ./mmsplice

RUN mkdir /root/.vep
RUN mkdir /root/.vep/Plugins/
RUN cp ./VEP_plugin/MMSplice.pm /root/.vep/Plugins/MMSplice.pm

RUN pip install -e .
WORKDIR ..
RUN mkdir outputs
RUN mmsplice

CMD ./vep
