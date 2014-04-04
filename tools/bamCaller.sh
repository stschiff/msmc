#!/usr/bin/env bash

SAMTOOLS=samtools
BCFTOOLS=bcftools
SAMTOOLS_OPT="-q 20 -Q 20 -C 50 -g"
BCFTOOLS_OPT="-cgI"
ZIP="bgzip -c"

usage() {
  echo "Usage: $0 [Options] <bam_file>
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

if [ -z $DEPTH ] || [ -z $CHR ]; then
  echo "need to specify depth" > /dev/stderr
  usage
fi

BASEDIR=$(dirname $0)
$SAMTOOLS mpileup $SAMTOOLS_OPT -r $CHR -f $REF $INPUT | $BCFTOOLS view $BCFTOOLS_OPT - | $BASEDIR/bamCallerHelper.py $DEPTH ${PREFIX}.mask.bed.gz | $ZIP > ${PREFIX}.vcf.gz

