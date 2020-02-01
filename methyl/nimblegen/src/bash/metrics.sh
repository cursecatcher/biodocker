#!/bin/bash

# rimuove il path e le estensioni dal nome del file
## funzione apposita per file con due estensioni !
function get_samplename {
  noext=${1##*/}    #leva path
  echo ${noext%%.*} #leva estensioni (2)
}

picard=/tools/picard/picard.jar
gatk=/tools/GATK/GenomeAnalysisTK.jar

alignment=$1
reference=$2
workspace=$3
metrics=$4
primary_target=$5
capture_target=$6

echo "INFO: indexing reference genome"
#un po' di variabili inutili
ref_directory=$(dirname $reference)
ref_name=$(basename $reference)
dict_ref=$ref_directory/${ref_name%.*}.dict

if [ -e $dict_ref ]; then
  rm $dict_ref
fi

java8 -jar $picard CreateSequenceDictionary REFERENCE=$reference OUTPUT=$dict_ref
samtools faidx $reference

echo "INFO: Indexing BAM files"
for bamfile in $(ls $alignment/*.bam); do
  samtools index $bamfile
done


cd $workspace

echo "INFO: basic mapping metrics"

for bamfile in $alignment/*.filtered.bam; do
  sample=`get_samplename $bamfile`
  metrics_file=$sample"_picard_alignment_metrics.txt"

  java8 -Xmx4g -Xms4g -jar $picard CollectAlignmentSummaryMetrics \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    INPUT=$bamfile \
    OUTPUT=$metrics_file \
    REFERENCE_SEQUENCE=$reference \
    VALIDATION_STRINGENCY=SILENT
done


echo "INFO: Hybrid Selection Metrics (todo)"

DESIGN_target_body=DESIGN_target_body.txt
DESIGN_bait_body=DESIGN_bait_body.txt

# Create a Picard Target Interval List Body
cat $primary_target | gawk '{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}' > $DESIGN_target_body
# Create a Picard Bait Interval List Body
cat $capture_target | gawk '{p$workspace_folder/rint $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}' > $DESIGN_bait_body

target_list=DESIGN_target_intervals.txt
bait_list=DESIGN_bait_intervals.txt

for bamfile in $alignment/*.clipped.bam; do
  sample=`get_samplename $bamfile`
  header_sample=$sample"_bam_header.txt"

  #create a picard interval list header
  samtools view -H $bamfile > $header_sample

  #concatenate header and target interval list to create a Picard Target Interval List
  cat $header_sample $DESIGN_target_body > $target_list
  cat $header_sample $DESIGN_bait_body > $bait_list

  hs_metrics=$sample"_picard_hs_metrics.txt"
  touch $hs_metrics #levare!

  java8 -Xmx4g -Xms4g -jar $picard CollectHsMetrics \
    BAIT_INTERVALS=$bait_list \
    TARGET_INTERVALS=$target_list \
    INPUT=$bamfile \
    OUTPUT=$hs_metrics \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    REFERENCE_SEQUENCE=$reference  \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=$workspace
done

echo "INFO: insert size metrics"

for bamfile in $alignment/*.filtered.bam; do
  sample=`get_samplename $bamfile`
  histogram=$sample"_picard_insert_size_plot.pdf"
  output=$sample"_insert_size_metrics.txt"

  java8 -Xmx4g -jar $picard  CollectInsertSizeMetrics \
    VALIDATION_STRINGENCY=LENIENT \
    HISTOGRAM_FILE=$histogram \
    INPUT=$bamfile \
    OUTPUT=$output
done


echo "INFO: Count on-target reads (todo)"
for bamfile in $alignment/*.filtered.bam; do
  sample=`get_samplename $bamfile`
  output_file=$sample"_on_target_reads_count.txt"

  #count on-target reads
  bedtools intersect -bed -abam $bamfile -b $primary_target > $output_file
done


echo "INFO: Calculate depth of coverage (todo)"

for bamfile in $alignment/*.clipped.bam; do
  sample=`get_samplename $bamfile`
  output_gatk=$sample"_gatk_primary_target_coverage"

 java8 -Xmx4g -Xms4g -jar $gatk -T DepthOfCoverage \
   -R $reference \
   -I $bamfile \
   -o $output_gatk \
   -L $primary_target \
   -ct 1 -ct 10 -ct 20
done


tar zcvf $metrics/picard_alignment_metrics.tar.gz *_picard_alignment_metrics.txt
tar zcvf $metrics/picard_hs_metrics.tar.gz *_picard_hs_metrics.txt
tar zcvf $metrics/picard_insert_size_plot.tar.gz *_picard_insert_size_plot.pdf
tar zcvf $metrics/insert_size_metrics.tar.gz *_insert_size_metrics.txt
tar zcvf $metrics/gatk_primary_target_coverage.tar.gz *_gatk_primary_target_coverage*
tar zcvf $metrics/on_target_reads_count.tar.gz *_on_target_reads_count.txt
rm -rf *
