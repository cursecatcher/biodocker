FROM ubuntu:16.04

ENV DEBIAN_FRONTEND noninteractive

# Install base system and devel tools
RUN \
  sed -i 's/# \(.*multiverse$\)/\1/g' /etc/apt/sources.list && \
  apt-get update && \
  apt-get install -y apt-utils && \
  apt-get -y upgrade && \
  apt-get install -y build-essential && \
  apt-get install -y software-properties-common && \
  apt-get install -y curl git htop man unzip vim wget && \
  apt-get install -y locales && \
  apt-get install -y libssl-dev libgit2-dev

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

#
# R external dependencies
#
RUN apt-get install -y openmpi-*
RUN apt-get install -y gfortran
RUN apt-get install -y xpdf
RUN apt-get install -y firefox
RUN apt-get install -y libreadline-dev
RUN apt-get install -y texlive-latex-base
RUN apt-get install -y texlive-fonts-recommended
RUN apt-get install -y texlive-fonts-extra
RUN apt-get install -y texlive-latex-extra
RUN apt-get install -y texinfo
RUN apt-get install -y xorg-dev
RUN apt-get install -y libbz2-dev
RUN apt-get install -y liblzma-dev
RUN apt-get install -y libpcre3-dev
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libtiff-dev
RUN apt-get install -y libjpeg-dev
RUN apt-get install -y libicu-dev
RUN apt-get install -y libcairo-dev
RUN apt-get install -y default-jdk
RUN apt-get install -y tk8.6-dev

# Compile and install an R in portable version
COPY files/R-3.5.1.tar.gz /tmp/R-3.5.1.tar.gz
WORKDIR /tmp
RUN tar zxf R-3.5.1.tar.gz
WORKDIR R-3.5.1
RUN ./configure --prefix=/opt/PortableR --enable-static=yes --enable-shared=yes --enable-R-shlib --with-tcltk

RUN make
RUN make install

WORKDIR /bin
RUN ln -s /opt/PortableR/bin/* .

WORKDIR /

#Setup CRAN mirror
COPY files/Rprofile.site /opt/PortableR/lib/R/etc

RUN /opt/PortableR/bin/R --version
RUN /opt/PortableR/bin/R -e "install.packages('littler')" && cd /bin && ln -s /opt/PortableR/lib/R/library/littler/bin/r .
RUN /opt/PortableR/bin/R -e "install.packages('devtools')"
RUN /opt/PortableR/bin/R -e "install.packages('BiocManager')"
RUN /opt/PortableR/bin/R -e "BiocManager::install('GenomicRanges')"
RUN apt-get install -y libxml2-dev
RUN /opt/PortableR/bin/R -e "BiocManager::install('BSgenome')"
RUN /opt/PortableR/bin/R -e "install.packages('stringi')"
RUN /opt/PortableR/bin/R -e "BiocManager::install('biomaRt')"
RUN /opt/PortableR/bin/R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')"

#Install bedtools
COPY files/bedtools-2.27.1.tar.gz /tmp/bedtools-2.27.1.tar.gz
WORKDIR /tmp
RUN tar zxf bedtools-2.27.1.tar.gz && \
    cd bedtools2 && \
    make && make install && make clean

RUN apt-get install -y gawk
#Install thread support for perl
RUN apt-get -y install cpanminus && cpanm threads

RUN apt-get -y install python3-pip && pip3 install pyliftover

# #Install samtools (not sure if it's required)
# COPY files/htslib-1.9.tar.bz2 /tmp/htslib-1.9.tar.bz2
# WORKDIR /tmp
# RUN tar xjf htslib-1.9.tar.bz2
# WORKDIR htslib-1.9
# RUN ./configure && make && make install
#
# COPY files/samtools-1.9.tar.bz2 /tmp/samtools-1.9.tar.bz2
# WORKDIR /tmp
# RUN tar xjf samtools-1.9.tar.bz2
# WORKDIR samtools-1.9
# RUN ./configure --without-curses && make && make install


#Install HashCirc 
COPY src/hashcirc/* /tmp/hashcirc/
RUN cd /tmp/hashcirc && make && make install 

#Install CIRI2, CIRI_AS and the relative python scripts 
COPY src/CIRI_AS/CIRI_AS_v1.2.pl src/CIRI_v2.0.6/CIRI2.pl src/docker4ciri/*.py  /ciri2/

#Install Circhunter 
COPY src/circhunter/* /scripts/circhunter/



WORKDIR /

# Create scratch directory
RUN mkdir /scratch; chmod 777 /scratch

# Create data directory
RUN mkdir /data; chmod 777 /data

# Tag the image with the actual repo commit
COPY BUILD_COMMIT /BUILD_COMMIT

# Dump all packages installed with it's version
RUN dpkg -l > /PKG_LIST

# Define default command.
CMD ["bash"]
