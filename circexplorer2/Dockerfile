FROM debian:latest 

#### https://www.jianshu.com/p/f584b6749b44

RUN apt-get -y update && apt-get -y upgrade 
RUN apt-get -y install gcc g++ make 
RUN apt-get -y install libbz2-dev zlib1g-dev libncurses5-dev libz-dev 
RUN apt-get -y install python-pip
RUN apt-get -y install liblzma-dev 
#install circexplorer2 
RUN pip install --upgrade setuptools && pip install pysam &&  pip install circexplorer2

COPY tools/tophat-2.1.1.Linux_x86_64.tar.gz \
     tools/cufflinks-2.2.1.Linux_x86_64.tar.gz \ 
     tools/bedtools-2.25.0.tar.gz \
     tools/2.7.0d.tar.gz \
     tools/bwa-0.7.17.tar.bz2 \
     tools/htslib-1.9.tar.bz2 \
     tools/segemehl-0.3.4.tar.gz \
     /tmp/
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
#install UCSC stuffs
COPY tools/bedGraphToBigWig \
     tools/bedToBigBed \ 
     tools/genePredToGtf \ 
     tools/gtfToGenePred \
     /usr/local/bin/
#install STAR 
WORKDIR /tmp
RUN tar xf 2.7.0d.tar.gz && cd STAR-2.7.0d/bin/Linux_x86_64_static && cp * /usr/local/bin/
#install bwa 
RUN tar xf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    cp bwa /usr/local/bin/
#install mapsplice 
COPY tools/MapSplice-v2.2.1.zip /tmp/
# RUN apt-get install -y unzip && \
#     unzip MapSplice-v2.2.1.zip && \
#     cd MapSplice-v2.2.1 && \
#     make && \
#     make install 
#install segemehl
RUN apt-get install -y libcurl4-openssl-dev libssl-dev pkg-config 
RUN tar xf htslib-1.9.tar.bz2 && cd htslib-1.9 && ./configure && make && make install 
RUN tar xf segemehl-0.3.4.tar.gz && \
    cd segemehl-0.3.4 && \
    export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig/:$PKG_CONFIG_PATH && \
    make

ENTRYPOINT [ "CIRCexplorer2" ]