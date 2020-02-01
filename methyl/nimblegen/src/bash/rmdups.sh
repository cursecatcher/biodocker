#!/bin/bash

# rimuove il path e le estensioni dal nome del file
## funzione apposita per file con una estensione !
function get_samplename {
  noext=${1##*/}    #leva path
  echo ${noext%.*} #leva estensione
}

picard=/tools/picard/picard.jar

alignments=$1
workspace=$2
metrics=$3

cd $workspace

for samfile in $(ls $alignments/*.sam); do
  sample=`get_samplename $samfile`

  # file_nopath=$(basename $samfile)
  # sample=${file_nopath%.*}

  bamfile=$sample.bam

  #Add read group + conversione sam2bam
  java8 -Xmx4g -Xms4g -jar $picard AddOrReplaceReadGroups \
      INPUT=$samfile  OUTPUT=$bamfile  \
      RGID=$sample RGLB=$sample RGPL=illumina RGSM=$sample RGPU=illumina

  #split by strand
  bamtools split -tag ZS -in $bamfile

  top_bam=$sample.top.bam
  bottom_bam=$sample.bottom.bam

  ### merge by strand - TOP
  top1=$sample".TAG_ZS_++.bam"
  top2=$sample".TAG_ZS_+-.bam"
  bamtools merge -in $top1 -in $top2 -out $top_bam

  ###merge by strand - BOTTOM
  bottom1=$sample".TAG_ZS_-+.bam"
  bottom2=$sample".TAG_ZS_--.bam"
  bamtools merge -in $bottom1 -in $bottom2 -out $bottom_bam

  rm $bamfile *.TAG_ZS_*

  top_sorted=$sample.top.sorted.bam
  bottom_sorted=$sample.bottom.sorted.bam

  samtools sort $top_bam -o $top_sorted
  samtools sort $bottom_bam -o $bottom_sorted

  rm $bottom_bam $top_bam

  ################################ remove duplicates
  top_nodups=$sample".top.rmdups.bam"
  bottom_nodups=$sample".bottom.rmdups.bam"

  metrics_top=$sample".top.rmdups_metrics.txt"
  metrics_bottom=$sample".bottom.rmdups_metrics.txt"

  bam_no_dups=$sample".rmdups.bam"

  java8 -Xmx4g -Xms4g -jar $picard MarkDuplicates \
      VALIDATION_STRINGENCY=LENIENT INPUT=$top_sorted  OUTPUT=$top_nodups \
      METRICS_FILE=$metrics_top  REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true

  java8 -Xmx4g -Xms4g -jar $picard MarkDuplicates \
      VALIDATION_STRINGENCY=LENIENT INPUT=$bottom_sorted  OUTPUT=$bottom_nodups \
      METRICS_FILE=$metrics_bottom  REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true

  bamtools merge -in $top_nodups -in $bottom_nodups -out $bam_no_dups
  ################################ end
  rm $top_sorted $bottom_sorted $top_nodups $bottom_nodups
  mv $bam_no_dups $bamfile
done

tar zcvf $metrics/rmdups_metrics.tar.gz *.rmdups_metrics.txt
rm *.txt *.bai
