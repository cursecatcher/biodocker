FROM debian:latest 

RUN apt update && apt upgrade -y 
RUN apt install -y build-essential cmake autoconf automake \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev \
    libssl-dev libncurses5-dev   \
    git wget python gawk curl   \
    r-base r-base-dev   \
    python3 python3-pip &&   pip3 install biopython 
RUN apt install -y nano tree

RUN mkdir /tools && mkdir home/data/ && mkdir /tools/picard && \
    wget -P /tools/picard https://github.com/broadinstitute/picard/releases/download/2.18.4/picard.jar 

RUN git clone https://github.com/samtools/htslib.git /tools/htslib && \ 
    git clone https://github.com/samtools/samtools.git /tools/samtools && \
    cd /tools/htslib && \
    autoheader && autoconf && ./configure && make && make install && \
    cd /tools/samtools && \
    autoheader && autoconf && ./configure && make && make install && \
    cd .. && rm -rf samtools htslib

RUN git clone https://github.com/arq5x/bedtools2.git /tools/bedtools2 && \
    cd /tools/bedtools2 && make && make install && cd .. && rm -rf bedtools2

RUN git clone git://github.com/pezmaster31/bamtools.git /tools/bamtools && \ 
    mkdir /tools/bamtools/build && cd /tools/bamtools/build && \
    cmake -DCMAKE_INSTALL_PREFIX=/tools/bamtools .. && make && make install &&\
    cp src/toolkit/bamtools /usr/local/bin/ \
    && cd /tools && rm -rf bamtools

RUN git clone https://github.com/statgen/bamUtil.git /tools/bamUtil && \ 
    cd /tools/bamUtil && make cloneLib && make && make install && \
    cd /tools && rm -rf bamUtil

RUN wget -P /tools http://lilab.research.bcm.edu/dldcc-web/lilab/yxi/bsmap/bsmap-2.90.tgz && \
    cd /tools && tar xf bsmap-2.90.tgz && cd bsmap-2.90 && \ 
    make && make install && cd /tools && rm -rf bsmap-2.90 bsmap-2.90.tgz

COPY pkgs/ /tmp/
WORKDIR /tmp 
RUN tar xf BisSNP.tar.gz && mv BisSNP /tools/ && \
    mkdir /tools/GATK && mv GenomeAnalysisTK.jar /tools/GATK  && \
    tar xf jre-7u80-linux-x64.tar.gz && \
    tar xf jre-8u172-linux-x64.tar.gz && \ 
    mkdir /usr/java && mv jre1.7.0_80 jre1.8.0_172 /usr/java/ && \
    ln -s /usr/java/jre1.7*/bin/java /usr/local/bin/java7 && \
    ln -s /usr/java/jre1.8*/bin/java /usr/local/bin/java8                                                                                                                            

COPY nimble.tar.gz /tmp/

RUN cd /tmp/ && tar zxvf nimble.tar.gz && \
    mv nimble.py /bin/ && mv bash/*.sh /bin/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

#ENTRYPOINT ["/bin/nimble.py"]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 