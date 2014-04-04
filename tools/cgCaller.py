#!/usr/bin/env python

import sys
import argparse
import gzip
import re
import string
import mask_generator

argparser = argparse.ArgumentParser()
argparser.add_argument("chr", help="Chromosome in the masterVar file")
argparser.add_argument("sample_id", help="Sample ID to put into the VCF file")
argparser.add_argument("out_prefix", help="Prefix for mask- and vcf file")
argparser.add_argument("input", help="Complete Genomics masterVarBeta file (uncompressed)")
argparser.add_argument("--max_pos", type=int, default=0)
args = argparser.parse_args()

mask_file = args.out_prefix + ".mask.bed.gz"
mask_generator = mask_generator.MaskGenerator(mask_file, args.chr)
vcf_file = gzip.open(args.out_prefix + ".vcf.gz", "w")

vcf_file.write("##fileformat=VCFv4.1\n")
vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(args.sample_id))

input_file = open(args.input, "r")
line_count = 0
chromosome_read = False
for line in input_file:
  if line[0] == '#' or line[0] == '>' or line == "\n":
    continue
  line_count += 1
  if line_count % 10000 == 0:
    sys.stderr.write("processing line {}\n".format(line_count))
    
  fields = string.split(string.strip(line))

  chrom = fields[2]
  begin = int(fields[3])
  end = int(fields[4])
  
  if chrom != args.chr:
    if chromosome_read:
      break
    else:
      continue

  if args.max_pos > 0 and end > args.max_pos:
      break
  chromosome_read = True
  zygosity = fields[5]
  var_type = fields[6]
  
  if var_type == "ref" and zygosity == "hom":
    for i in xrange(begin + 1, end + 1):
      mask_generator.addCalledPosition(i)

  if var_type == "snp":
    if zygosity in ["hom", "het-ref", "het-alt"]:
      assert end - begin == 1
      allele_ref = fields[7]
      allele_1 = fields[8]
      allele_2 = fields[9]
      allele1_qual = fields[14]
      allele2_qual = fields[15]
      if allele1_qual == "VQHIGH" and allele2_qual == "VQHIGH":
        mask_generator.addCalledPosition(begin + 1)
        allele_indices = []
        alt_alleles = []
        if allele_1 != allele_ref:
          alt_alleles.append(allele_1)
          allele_indices.append(1)
        if allele_2 == allele_ref:
          allele_indices.append(0)
        elif allele_2 == allele_1:
          allele_indices.append(1)
        else:
          alt_alleles.append(allele_2)
          allele_indices.append(2)
        vcf_file.write("{chrom}\t{pos}\t.\t{ref_a}\t{alt_a}\t.\tPASS\t.\tGT\t{gen1}/{gen2}\n".format(chrom=args.chr, 
                       pos=begin+1, ref_a=allele_ref, alt_a=",".join(alt_alleles), gen1=allele_indices[0], gen2=allele_indices[1]))
