FROM debian:latest 

#### https://www.jianshu.com/p/f584b6749b44

RUN apt-get -y update && apt-get -y upgrade && apt-get -y install \
    gcc g++ make \
    libbz2-dev zlib1g-dev libncurses5-dev libz-dev liblzma-dev \
    python-pip
#install circexplorer2 
RUN pip install --upgrade setuptools && pip install pysam &&  pip install circexplorer2

COPY tools/tophat-2.1.1.Linux_x86_64.tar.gz \
     tools/cufflinks-2.2.1.Linux_x86_64.tar.gz \ 
     tools/bedtools-2.25.0.tar.gz \
     tools/htslib-1.9.tar.bz2 \
     /tmp/
#     tools/2.7.0d.tar.gz \
#     tools/bwa-0.7.17.tar.bz2 \
#     tools/segemehl-0.3.4.tar.gz \
#install UCSC stuffs
COPY tools/bedGraphToBigWig \
     tools/bedToBigBed \ 
     tools/genePredToGtf \ 
     tools/gtfToGenePred \
     /usr/local/bin/


WORKDIR /tmp
#install top-hat & top-hat fusion 
RUN tar xf tophat-2.1.1.Linux_x86_64.tar.gz && \
    cd tophat-2.1.1.Linux_x86_64 && \
    cp b* c* f* g* j* long_spanning_reads map2gtf prep_reads sam* segment_juncs sra_to_solid tophat* /usr/local/bin/
#install cufflinks
RUN tar xf cufflinks-2.2.1.Linux_x86_64.tar.gz && \
    cd cufflinks-2.2.1.Linux_x86_64 && \
    cp * /usr/local/bin 
#install bedtools 
RUN tar xf bedtools-2.25.0.tar.gz && cd bedtools2 && make && make install 

#install STAR 
#WORKDIR /tmp
#RUN tar xf 2.7.0d.tar.gz && cd STAR-2.7.0d/bin/Linux_x86_64_static && cp * /usr/local/bin/
#install bwa 
#RUN tar xf bwa-0.7.17.tar.bz2 && \
#    cd bwa-0.7.17 && \
#    make && \
#    cp bwa /usr/local/bin/
#install mapsplice (anche no)
#COPY tools/MapSplice-v2.2.1.zip /tmp/
# RUN apt-get install -y unzip && \
#     unzip MapSplice-v2.2.1.zip && \
#     cd MapSplice-v2.2.1 && \
#     make && \
#     make install 
#install segemehl
#RUN apt-get install -y libcurl4-openssl-dev libssl-dev pkg-config 
#RUN tar xf htslib-1.9.tar.bz2 && cd htslib-1.9 && ./configure && make && make install 
#RUN tar xf segemehl-0.3.4.tar.gz && \
#    cd segemehl-0.3.4 && \
#    export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig/:$PKG_CONFIG_PATH && \
#    make

#COPY tools/bowtie-1.2.3-linux-x86_64.zip tools/bowtie2-2.3.5.1-linux-x86_64.zip /tmp/
#RUN apt install -y unzip 
#RUN unzip bowtie2-2.3.5.1-linux-x86_64.zip && mv bowtie2-2.3.5.1-linux-x86_64/bowtie2* /usr/local/bin/ && rm -r bowtie2-2.3.5.1-linux-x86_64 
#RUN unzip bowtie-1.2.3-linux-x86_64.zip && mv bowtie-1.2.3-linux-x86_64/bowtie* /usr/local/bin && rm -r bowtie-1.2.3-linux-x86_64

COPY src/circexplorer2.py /bin/

#RUN tar xf /bin/bowtie-build.tar.gz 
# && rm -f /bin/bowtie-build.tar.gz 

ENTRYPOINT [ "/bin/circexplorer2.py" ]