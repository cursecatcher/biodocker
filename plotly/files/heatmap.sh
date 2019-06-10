#!/bin/bash

#filename="${1}"
color_palette="${1}"
separator="${2}"
status="${3}"
#lower_bound="${4}"
#upper_bound="${5}"
logfile="/data/out/heatmaply.log"

scratch="/scratch"

##### Creating a scratch folder for the current task
scratch_id=$(date '+%d-%m-%Y-%H-%M-%S-%N')
my_scratch=${scratch}/${scratch_id}
echo "Creating scratch folder ${my_scratch}" &>> ${logfile}
mkdir ${my_scratch}


echo "Preprocessing count table: filtering for selected genes" &>> ${logfile}
python3 /4seq/preprocessing.py ${my_scratch} ${separator} ${status} &>> ${logfile}
exitvalue=$?

if [ $exitvalue -ne 0 ]; then
    echo "Count table preprocessing failed. Abort"
    exit $exitvalue
fi

echo "Creating heatmap"
Rscript /4seq/heatmap.R ${my_scratch} ${color_palette} &>> ${logfile}
#Rscript /4seq/heatmap.R ${my_scratch} ${filename} ${lower_bound} ${upper_bound} &>> ${logfile}
exitvalue=$?

echo "Execution completed. Deleting temporary files and exiting" &>> ${logfile}
rm -rf ${my_scratch}

exit $exitvalue
