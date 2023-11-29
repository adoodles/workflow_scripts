FROM ubuntu:20.04

RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y python python-dev \
    python3 python3-dev python3-pip apt-transport-https openjdk-8-jre wget \
    zip git-all zlib1g-dev libbz2-dev liblzma-dev libjpeg-dev

RUN pip3 install cloudpickle

RUN pip3 install anadama2 

RUN apt-get install -y bowtie2

RUN pip3 install Cython

RUN pip3 install numpy

RUN pip3 install biom-format

RUN git clone https://github.com/SegataLab/cmseq.git && \
    cd cmseq && python3 setup.py install && cd ../ && rm -r cmseq

RUN pip3 install metaphlan==4.0.6

WORKDIR /tmp