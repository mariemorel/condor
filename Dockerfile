
FROM ubuntu:22.04

ENV SINGULARITY_VERSION="v3.10.4"
ENV SINGULARITY_DEB="singularity-ce_3.10.4-jammy_amd64.deb"

COPY . /usr/local/bin/
RUN apt-get update --fix-missing && \
    apt-get install -y wget && \
    apt-get install -y libglib2.0-0 cryptsetup-bin squashfs-tools runc ca-certificates-java openjdk-18-jre curl && \
    useradd -ms /bin/bash condor && \
    cd /usr/local/ && \
    wget https://github.com/sylabs/singularity/releases/download/${SINGULARITY_VERSION}/${SINGULARITY_DEB} && \
    dpkg -i ${SINGULARITY_DEB} && \
    rm ${SINGULARITY_DEB} && \
    mkdir /usr/local/singularity && \
    cd /usr/local/bin && curl -s https://get.nextflow.io | bash && ./nextflow -h && \
    chmod 777 /usr/local/singularity && \
    chmod 777 /usr/local/bin/nextflow && \
    chmod 777 /usr/local/bin/run_condor

USER condor

ENV NXF_TEMP=/home/condor/tmp/
RUN mkdir /home/condor/tmp/ && nextflow -h && \
    singularity pull --name /usr/local/singularity/evolbioinfo-python-evol-v3.7.3condor.img  docker://evolbioinfo/python-evol:v3.7.3condor && \
    singularity pull --name /usr/local/singularity/evolbioinfo-gotree:v0.3.0b.img docker://evolbioinfo/gotree:v0.3.0b && \
    singularity pull --name /usr/local/singularity/evolbioinfo-goalign-v0.3.3c.img docker://evolbioinfo/goalign:v0.3.3c && \
    singularity pull --name /usr/local/singularity/evolbioinfo-iqtree-v2.1.3.img docker://evolbioinfo/iqtree:v2.1.3 && \
    singularity pull --name /usr/local/singularity/evolbioinfo-iqtree-v1.6.8.img docker://evolbioinfo/iqtree:v1.6.8 && \
    singularity pull --name /usr/local/singularity/evolbioinfo-pastml-v1.9.33.img docker://evolbioinfo/pastml:v1.9.33 && \
    singularity pull --name /usr/local/singularity/evolbioinfo-bayestraits-v3.0.1.img docker://evolbioinfo/bayestraits:v3.0.1

ENTRYPOINT ["/usr/local/bin/run_condor"]
