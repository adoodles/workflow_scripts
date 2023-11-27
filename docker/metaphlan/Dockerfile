FROM ubuntu:20.04

RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y python python-dev \
    python3 python3-dev python3-pip apt-transport-https openjdk-8-jre wget \
    zip git-all zlib1g-dev libbz2-dev liblzma-dev libjpeg-dev

RUN pip3 install cloudpickle

RUN pip3 install anadama2 

RUN mkdir -p $HOME/tools/bowtie2/ && cd $HOME/tools/bowtie2/ && \
    wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip/download \
    && unzip download && export PATH=$HOME/tools/bowtie2/bowtie2-2.4.2-sra-linux-x86_64:$PATH && cd /

RUN pip3 install Cython

RUN pip3 install numpy

RUN pip3 install biom-format

RUN git clone https://github.com/SegataLab/cmseq.git && \
    cd cmseq && python3 setup.py install && cd ../ && rm -r cmseq

RUN pip3 install metaphlan==4.0.6

WORKDIR /tmp