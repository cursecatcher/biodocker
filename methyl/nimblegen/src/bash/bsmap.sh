#!/bin/bash

sample=$1
sample_r1=$2
sample_r2=$3
genome=$4
alignment=$5
#workspace=$6

# picard=/tools/picard/picard.jar

samfile=$alignment/$sample.sam
bamfile=$alignment/$sample.bam

#Alignment
bsmap -r 0 -s 16 -n 1 -a $sample_r1 -b $sample_r2 -d $genome -o $samfile #" #  2> $sample_r1.log"

exit 0
