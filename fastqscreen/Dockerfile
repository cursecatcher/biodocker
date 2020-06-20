FROM stevetsa/fastqscreen:0.12.0 

WORKDIR /tmp
COPY pkgs/ /tmp/

RUN tar xf Bismark-0.22.3.tar.gz && \
    cd Bismark-0.22.3 && \
    cp *bismark* coverage2cytosine filter_non_conversion methylation_consistency NOMe_filtering /usr/local/bin/ 

RUN tar xf v0.14.1.tar.gz && \
    cd FastQ-Screen-0.14.1 && \
    cp fastq_screen /usr/local/bin

RUN cp fastq_screen.conf /usr/local/bin/


WORKDIR /