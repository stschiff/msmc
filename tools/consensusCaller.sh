#!/usr/bin/env bash

usage() {
  echo "Usage: $0 [Options] <input_file>
input_file: can be either a bam file or a vcf file in Complete Genomics format (ending *.vcf.gz or *.vcf)
Options:
-r <ref_file>: the reference file
-c <chr>: the chromosome
-d <mean_depth>: the mean read depth
-o <out_prefix>" > /dev/stderr
  exit 1
}

MODE="BAM"
while getopts ":r:c:d:o:" opt; do
  case $opt in
    r)
    REF=$OPTARG
    ;;
    c)
    CHR=$OPTARG
    ;;
    d)
    DEPTH=$OPTARG
    ;;
    o)
    PREFIX=$OPTARG
    ;;
    *)
    usage
    ;;
  esac
done
shift $((OPTIND-1))

INPUT=$1
TYPE=$(echo $INPUT | awk '{l=length($0); print substr($0, l - 2,3)}')

if [ ! -f $REF ] || [ -z $REF ]; then
  echo "did not find reference file" > /dev/stderr
  usage
fi

if [ ! -f $INPUT ] || [ -z $INPUT ]; then
  echo "did not find bam file" > /dev/stderr
  usage
fi

if [ -z $PREFIX ]; then
  echo "need to specify output prefix" > /dev/stderr
  usage
fi

BASEDIR=$(dirname $0)

case $TYPE in 
  bam)
  if [ -z $DEPTH ] || [ -z $CHR ]; then
    echo "need to specify depth" > /dev/stderr
    usage
  fi
  samtools mpileup -q 20 -Q 20 -C 50 -g -r $CHR -f $REF $INPUT | bcftools view -cgI - | $BASEDIR/vcfParser.py $DEPTH ${PREFIX}.mask.bed.gz "BAM" $REF $CHR | bgzip -c > ${PREFIX}.vcf.gz
  ;;
  vcf)
  DEPTH=10 # use dummy depth of 10, ignored anyway
  awk -v c=$CHR '$0 ~ /^#/ || $1 == c' $INPUT | $BASEDIR/vcfParser.py $DEPTH ${PREFIX}.mask.bed.gz "CG" $REF $CHR | bgzip -c > ${PREFIX}.vcf.gz
  ;;
  .gz)
  DEPTH=10 # use dummy depth of 10, ignored anyway
  cat $INPUT | gzip -d | awk -v c=$CHR '$0 ~ /^#/ || $1 == c' | $BASEDIR/vcfParser.py $DEPTH ${PREFIX}.mask.bed.gz "CG" $REF $CHR | bgzip -c > ${PREFIX}.vcf.gz
  ;;
  *)
  echo "unknown file type, need BAM or Complete Genomics VCF" > /dev/stderr
  ;;
esac

