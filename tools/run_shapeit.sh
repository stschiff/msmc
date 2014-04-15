#!/usr/bin/env bash

VCF=$1
TMP_DIR=$2
CHR=$3

echo "VCF: $VCF" > /dev/stderr
echo "TMP_DIR: $TMP_DIR" > /dev/stderr
echo "CHR: $CHR" > /dev/stderr

mkdir -p $TMP_DIR

REF_DIR=$lhome/1KG_impute/2013-09-16
REF_PANEL="$REF_DIR/ALL.chr$CHR.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz $REF_DIR/ALL.chr$CHR.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz $REF_DIR/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample"
MAP=$REF_DIR/genetic_map_chr${CHR}_combined_b37.txt

NAME=$(basename $VCF .vcf.gz)
DIR=$(dirname $VCF)
TMP_VCF=$TMP_DIR/$NAME.noMultiAllelicSites.vcf.gz
PHASED_VCF=$TMP_DIR/$NAME.onlyPhased.vcf
SHAPEIT_OUT="$TMP_DIR/$NAME.phased.haps.gz $TMP_DIR/$NAME.phased.samples"
MERGED_VCF=$DIR/$NAME.phased.vcf.gz

bcftools-exp-rc view -M 2 -O z $VCF > $TMP_VCF

# This step is necessary to make a list of sites to exclude from the main run:
shapeit -check -V $TMP_VCF -M $MAP --input-ref $REF_PANEL --output-log $TMP_DIR/$NAME.alignments

# Main run
shapeit -V $TMP_VCF -M $MAP --input-ref $REF_PANEL -O $SHAPEIT_OUT --exclude-snp $TMP_DIR/$NAME.alignments.snp.strand.exclude --no-mcmc --output-log $TMP_DIR/$NAME.main

# Convert back to phased vcf, containing only the phased sites that shapeit used
shapeit -convert --input-haps $SHAPEIT_OUT --output-vcf $PHASED_VCF --output-log $TMP_DIR/$NAME.convert

# zipping and indexing
bcftools-exp-rc view -O z $PHASED_VCF > $PHASED_VCF.gz
bcftools-exp-rc index -f $PHASED_VCF.gz
bcftools-exp-rc index -f $VCF

# merging phased and unphased vcfs, keeping all unphased sites from the original vcf, but replacing the phased ones.
bcftools-exp-rc merge $VCF $PHASED_VCF.gz | awk '
  BEGIN {OFS="\t"}
  $0 ~ /^##/ {print}
  $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}
  $0 !~ /^#/ {
    if(substr($11, 1, 3) != "./.")
      $10 = $11
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
  }' | bcftools-exp-rc view -O z > $MERGED_VCF
