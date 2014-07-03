#!/usr/bin/env bash

# BAM Calling
# REF=$lhome/ref/hs37d5.fa
# REGION="20:10000000-10100000"
# COV=30
# BAM=/lustre/scratch113/projects/t19-ethiopia/vrpipe/results/6/c/5/1/279885/4_bam_mark_duplicates/pe.1.markdup.bam
#
# samtools-exp-rc mpileup -q 20 -Q 20 -C 50 -u -f $REF -r $REGION $BAM | bcftools-exp-rc call -c -V indels | ~/Projects/msmc/tools/bamCaller.py $COV test.bed.gz | gzip -c > test.vcf.gz

