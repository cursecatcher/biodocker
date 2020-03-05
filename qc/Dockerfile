FROM piegsaj/oracle-jre:1.8.0_172-b11

RUN apt update && apt upgrade -y 
RUN apt install -y python3 python3-pip && pip3 install multiqc 
COPY pkgs/FastQC /opt/FastQC
WORKDIR /opt/FastQC
RUN chmod +x fastqc && \
    ln -s /opt/FastQC/fastqc /usr/bin/fastqc 
COPY start_qc.sh /usr/bin/
WORKDIR /

ENTRYPOINT [ "start_qc.sh" ]