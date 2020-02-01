#!/bin/bash

# rimuove il path e le estensioni dal nome del file
## funzione apposita per file con due estensioni !
function get_samplename {
  noext=${1##*/}    #leva path
  echo ${noext%%.*} #leva estensioni (2)
}

bissnp=/tools/BisSNP/BisSNP.jar
bissnp_utils=/tools/BisSNP/utils/

alignment=$1
genome=$2
workspace=$3
results=$4
dbsnp=$5
capture_target=$6

reference_index=$genome.fai

cd $workspace

for bamfile in $alignment/*.clipped.bam; do
  sample=`get_samplename $bamfile`

  recal_file=$sample.recalFile_before.csv
  recal_bam=$sample.recal.bam

  ### base quality recalibration
  java7 -Xmx10g -jar $bissnp \
    -R $genome \
    -I $bamfile \
    -T BisulfiteCountCovariates \
    -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate \
    -recalFile $recal_file \
    -knownSites $dbsnp \
    -nt 8 #num processors

    java7 -Xmx10g -jar $bissnp \
      -R $genome \
      -I $bamfile \
      -o $recal_bam \
      -T BisulfiteTableRecalibration \
      -recalFile $recal_file \
      -maxQ 40

    ### combined SNP and methylation calling
    cpg_raw_vcf=$sample.cpg.raw.vcf
    snp_raw_vcf=$sample.snp.raw.vcf

    java7 -Xmx10g -jar $bissnp \
      -R $genome \
      -I $recal_bam \
      -T BisulfiteGenotyper \
      -D $dbsnp \
      -vfn1 $cpg_raw_vcf -vfn2 $snp_raw_vcf \
      -L $capture_target \
      -stand_call_conf 20 -stand_emit_conf 0 -mmq 30 -mbq 0 \
      -nt 8 #numero processori

    ###sort vcf files
    cpg_sorted_vcf=$sample.cpg.raw.sorted.vcf
    snp_sorted_vcf=$sample.snp.raw.sorted.vcf

    ##scriptino perl x 2
    perl $bissnp_utils/sortByRefAndCor.pl --k 1 --c 2 $snp_raw_vcf $reference_index > $snp_sorted_vcf
    perl $bissnp_utils/sortByRefAndCor.pl --k 1 --c 2 $cpg_raw_vcf $reference_index > $cpg_sorted_vcf


    #####Filter SNP/methylation calls
    #output filenames
    snp_def_vcf=$sample.snp.filtered.vcf
    cpg_def_vcf=$sample.cpg.filtered.vcf

    output_snp=$sample.snp.filter.summary.txt
    output_cpg=$sample.cpg.filter.summary.txt

    java7 -Xmx10g -jar $bissnp \
      -R $genome \
      -T VCFpostprocess \
      -oldVcf $snp_sorted_vcf \
      -newVcf $snp_def_vcf \
      -snpVcf $snp_sorted_vcf \
      -o $output_snp

    java7 -Xmx10g -jar $bissnp \
      -R $genome \
      -T VCFpostprocess \
      -oldVcf $cpg_sorted_vcf \
      -newVcf $cpg_def_vcf \
      -snpVcf $snp_sorted_vcf \
      -o $output_cpg

      #VCF2bed conversion
      ####perl script
      perl $bissnp_utils/vcf2bed6plus2.strand.pl $snp_def_vcf
      perl $bissnp_utils/vcf2bed6plus2.strand.pl $cpg_def_vcf
done

tar zcvf $results/cpg_bed.tar.gz *.cpg.*.bed
tar zcvf $results/snp_bed.tar.gz *.snp.*.bed
tar zcvf $results/cpg_filter_summary.tar.gz *.cpg.filter.summary.txt
tar zcvf $results/snp_filter_summary.tar.gz *.snp.filter.summary.txt
rm -f *

exit 0
