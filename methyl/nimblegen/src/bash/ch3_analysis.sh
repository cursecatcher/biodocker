#!/bin/bash

# rimuove il path e le estensioni dal nome del file
## funzione apposita per file con due estensioni !
function get_samplename {
  noext=${1##*/}    #leva path
  echo ${noext%%.*} #leva estensioni (2)
}

samtools_path=/tools/BSMAP/samtools/

alignment=$1
genome=$2
workspace=$3
control_chr=$4
results=$5

cd $workspace

for bamfile in $alignment/*.clipped.bam; do
  sample=`get_samplename $bamfile`
  output=$sample.methylation_results.bed

  methratio.py -d $genome -s $samtools_path  -m 1 -z -i skip \
    -o $output $bamfile
done

echo "INFO: bisulfite conversion efficiency estimation"

for bamfile in $alignment/*.clipped.bam; do
  sample=`get_samplename $bamfile`
  output=$sample.efficiency_results.bed

  methratio.py -d $genome -s $samtools_path -m 1 -z -i skip \
    -o $output -c $control_chr $bamfile
done

tar zcvf $results/methylation_results.tar.gz *.methylation_results.bed
tar zcvf $results/efficiency_results.tar.gz *.efficiency_results.bed
rm -f * 

exit 0
