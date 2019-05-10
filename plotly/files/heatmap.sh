#!/bin/bash 

filename="${1}"
separator="${2}"
status="${3}"
lower_bound="${4}"
upper_bound="${5}"

scratch="/scratch"

##### Creating a scratch folder for the current task
scratch_id=$(date '+%d-%m-%Y-%H-%M-%S-%N')
my_scratch=${scratch}/${scratch_id}
echo "Creating scratch folder ${my_scratch}"
mkdir ${my_scratch}


echo "Preprocessing count table: filtering for selected genes"
python3 /4seq/preprocessing.py ${my_scratch} ${separator} ${status}
echo "Creating heatmap"
Rscript /4seq/heatmap.R ${my_scratch} ${filename} ${lower_bound} ${upper_bound}
chmod 666 /data/out/*

echo "Execution completed. Deleting temporary files and exiting"
rm -rf ${my_scratch}
