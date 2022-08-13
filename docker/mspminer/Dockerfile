FROM ubuntu:20.04

RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y python python-dev python3 python3-dev python3-pip wget unzip

RUN mkdir -p Downloads; \
    wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh -P /Downloads; \
    chmod +rwx /Downloads/Anaconda3-2022.05-Linux-x86_64.sh;

RUN echo "a7c0afe862f6ea19a596801fc138bde0463abcbce1b753e8d5c474b506a2db2d /Downloads/Anaconda3-2022.05-Linux-x86_64.sh" \
    | sha256sum -c

RUN yes "yes" | bash /Downloads/Anaconda3-2022.05-Linux-x86_64.sh

RUN . ~/.bashrc && conda install -c bioconda megahit

RUN wget https://www.enterome.com/wp-content/uploads/2022/01/mspminer_bin.zip -P /Downloads; \
    mkdir -p msp; \
    unzip /Downloads/mspminer_bin.zip -d /msp;

RUN DEBIAN_FRONTEND="noninteractive" apt-get install -y bwa samtools prodigal cd-hit
