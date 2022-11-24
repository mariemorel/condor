
FROM ubuntu:22.04

ENV SINGULARITY_VERSION="v3.10.4"
ENV SINGULARITY_DEB="singularity-ce_3.10.4-jammy_amd64.deb"

COPY . /usr/local/bin/
RUN apt-get update --fix-missing && \
    apt-get install -y wget && \
    apt-get install -y libglib2.0-0 cryptsetup-bin squashfs-tools runc ca-certificates-java openjdk-18-jre curl && \
    cd /usr/local/ && \
    wget https://github.com/sylabs/singularity/releases/download/${SINGULARITY_VERSION}/${SINGULARITY_DEB} && \
    dpkg -i ${SINGULARITY_DEB} && \
    rm ${SINGULARITY_DEB} && \
    mkdir /usr/local/singularity && \
    cd /usr/local/bin && curl -s https://get.nextflow.io | bash && ./nextflow -h

ENTRYPOINT ["/usr/local/bin/run_condor"]
