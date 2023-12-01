FROM ubuntu:20.04

RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y python python-dev \
    python3 python3-dev python3-pip build-essential cmake \
    openjdk-8-jre wget zip ncbi-blast+ git-all zlib1g-dev libffi-dev libbz2-dev liblzma-dev libjpeg-dev \
    racon 

RUN apt-get install -y spades

RUN git clone https://github.com/rrwick/Unicycler.git; \
    cd Unicycler; \
    python3 setup.py install; cd / \

WORKDIR /tmp