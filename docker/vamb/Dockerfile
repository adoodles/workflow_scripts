FROM ubuntu:20.04

RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y python python-dev python3 python3-dev python3-pip wget

RUN apt-get install -y libbz2-dev liblzma-dev

RUN pip3 install Cython

RUN pip3 install https://github.com/RasmussenLab/vamb/archive/v3.0.3.zip

RUN DEBIAN_FRONTEND="noninteractive" apt-get install -y spades minimap2 samtools
