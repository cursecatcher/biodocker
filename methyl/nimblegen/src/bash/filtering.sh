#!/bin/bash

workspace=$1
alignments=$2

cd $workspace

for bamfile in $(ls *.bam); do
  sample=${bamfile%.*}
  filtered_file=$sample".filtered.bam"
  clipped_file=$sample".clipped.bam"

  bamtools filter -isMapped true -isPaired true -isProperPair true -forceCompression \
    -in $bamfile -out $filtered_file

  bam clipOverlap --stats --in $filtered_file --out $clipped_file

  rm $bamfile
done

mv *.filtered.bam *.clipped.bam $alignments/
