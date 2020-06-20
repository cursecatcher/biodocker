FROM piegsaj/oracle-jre:1.8.0_172-b11

RUN apt update && apt upgrade -y 
RUN apt install -y \
    pigz \
    python3 \
    python3-pip && \
    pip3 install multiqc cutadapt

COPY pkgs/FastQC /opt/FastQC
RUN chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/bin/fastqc 
COPY src/ /usr/bin/
WORKDIR /

#ENTRYPOINT [ "start_qc.sh" ]
ENTRYPOINT [ "qc.py" ]
CMD [ "-h" ]