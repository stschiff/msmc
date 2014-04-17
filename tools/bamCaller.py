#!/usr/bin/env python3

# Use as: samtools mpileup -q 20 -Q 20 -C 50 -g -r <chr> -f <ref_file.fa> <input.bam> | bcftools view -cgI - | bamCaller.py <depth> <out_mask> > <out_vcf>
# or with new samtools:
# samtools-exp-rc mpileup -q 20 -Q 20 -C 50 -g -r <chr> -f <ref_file.fa> <input.bam> | bcftools-exp-rc call -c -S indels | bamCaller.py <depth> <out_mask> [Options...] | bcftools-exp-rc view -O z > <out_vcf.gz>

import sys
import re
import gzip
import utils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("depth", type=float)
parser.add_argument("mask")
parser.add_argument("--minMapQ", type=float, default=20)
parser.add_argument("--minConsQ", type=float, default=20)
parser.add_argument("--legend_file", help="Impute2 reference panel legend file, can be gzipped or not")
args = parser.parse_args()

sites_parser = None
if args.legend_file is not None:
    sites_parser = utils.LegendParser(args.legend_file)

minDepth = args.depth / 2.0
maxDepth = args.depth * 2.0

lastPos = 0
line_cnt = 0
chr_ = ""
for line in sys.stdin:
    if line[0] == '#':
        print(line, end="")
        continue
    fields = line.strip().split('\t')
    if chr_ == "":
        chr_ = fields[0]
        mask = utils.MaskGenerator(args.mask, chr_)
    else:
        assert fields[0] == chr_, "found multiple chromosomes"
    pos = int(fields[1])
    refAllele = fields[3]
    altAllele = fields[4]
    info = fields[7]
    genotypes = fields[9]
    if line_cnt % 10000 == 0:
        print("parsing position {}".format(pos), file=sys.stderr)
    line_cnt += 1

    if sites_parser is not None:
        while not sites_parser.end and sites_parser.pos < pos:
            sites_parser.tick()
    
    if re.match("^[ACTGactg]$", refAllele) and re.match("^[ACTGactg\.]$", altAllele):
        dp_match = re.search("DP=(\d+)", info)
        mq_match = re.search("MQ=(\d+)", info)
        fq_match = re.search("FQ=([\d-]+)", info)
        if not (dp_match and mq_match and fq_match):
            continue
        dp = int(dp_match.group(1))
        mq = int(mq_match.group(1))
        fq = int(fq_match.group(1))
        if dp >= minDepth and dp <= maxDepth and mq >= args.minMapQ and abs(fq) >= args.minConsQ:
            mask.addCalledPosition(pos)
            if altAllele != "." and re.match("^[01][/|][01]", genotypes):
                print(line, end="")
            elif sites_parser is not None and sites_parser.pos == pos:
                assert refAllele == sites_parser.ref_a
                fields[4] = sites_parser.alt_a
                print("\t".join(fields))

