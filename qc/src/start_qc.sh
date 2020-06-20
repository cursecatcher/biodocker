#!/bin/bash 

#to fix encoding stuff of multiqc using python3 
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

ifolder="/data"
ofolder="${ifolder}/quality_control"

#todo - check $ofolder existence?
mkdir -p ${ofolder}

echo "Starting FastQC with the following options: $*"
fastqc -o ${ofolder} $* ${ifolder}/*.gz ${ifolder}/*.fastq 

if [ $? -eq 0 ]; then
    cd ${ofolder}
    echo "Starting MultiQC..."
    multiqc . 
    chmod 777 ${ofolder} ${ofolder}/*
else
    echo "FastQC terminated with errors. Check input parameters and try again."
    rm -rf ${ofolder}
fi 