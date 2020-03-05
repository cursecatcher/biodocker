#!/bin/bash 

#to fix a strange behaviour of multiqc using python3 
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

ifolder="/data"
ofolder="${ifolder}/quality_control"
nthreads="${1}"

if [ -z ${nthreads+x} ]; then
    nthreads=1
fi

mkdir -p ${ofolder}
cd ${ofolder}

echo "Starting FastQC..."
fastqc --noextract -o ${ofolder} -t ${nthreads} -f fastq ${ifolder}/*.gz ${ifolder}/*.fastq &> ${ofolder}/fastqc.log 


echo "Starting MultiQC..."
multiqc . 

chmod 777 ${ofolder} ${ofolder}/*