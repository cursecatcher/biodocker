FROM comics/bowtie2:2.3.4.1

RUN yum install -y \
    make \
    lbzip2 \
    libcurl-devel \
    zlib-devel \
    bzip2-devel \
    xz-devel 
    
#RUN apt install -y make zlib-dev bzip2-dev xz-dev
COPY pkgs/ /tmp/
WORKDIR /tmp
RUN tar xf samtools-1.10.tar.bz2 && \
    cd samtools-1.10 && \
    ./configure --without-curses && \
    make && make install && \ 
    cd .. && rm -rf samtools-1.10 

#COPY Bismark-0.22.3.tar.gz /tmp/
RUN tar xf Bismark-0.22.3.tar.gz && \
    cd Bismark-0.22.3 && \
    cp *bismark* coverage2cytosine filter_non_conversion methylation_consistency NOMe_filtering /usr/local/bin/ 

#COPY bsmap-2.90.tgz /tmp/
RUN tar xf bsmap-2.90.tgz && \ 
    cd bsmap-2.90 && \
    make && \
    make install && \ 
    cd ..  && rm -rf bsmap-2.90
WORKDIR /

RUN yum install -y python3 
COPY src/ /usr/bin/

ENTRYPOINT [ "bsmark.py" ]
CMD [ "-h" ]