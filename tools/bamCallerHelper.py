#!/usr/bin/env python

import sys
import re
import gzip
import mask_generator

minMapQ = 20;
minConsQ = 20;

mean_depth = float(sys.argv[1])
mask_filename = sys.argv[2]

minDepth = mean_depth / 2.0
maxDepth = mean_depth * 2.0

lastPos = 0
line_cnt = 0
chr_ = ""
for line in sys.stdin:
  if line[0] == '#':
    print line,
    continue
  fields = line.strip().split('\t')
  if chr_ == "":
    chr_ = fields[0]
    mask = mask_generator.MaskGenerator(mask_filename, chr_)
  else:
    assert fields[0] == chr_, "found multiple chromosomes"
  pos = int(fields[1])
  refAllele = fields[3]
  altAllele = fields[4]
  info = fields[7]
  genotypes = fields[9]
  if line_cnt % 10000 == 0:
    sys.stderr.write("parsing position {}\n".format(pos))
  line_cnt += 1

  if re.match("^[ACTGactg]$", refAllele) and re.match("^[ACTGactg\.]$", altAllele):
    dp_match = re.search("DP=(\d+)", info)
    mq_match = re.search("MQ=(\d+)", info)
    fq_match = re.search("FQ=([\d-]+)", info)
    if not (dp_match and mq_match and fq_match):
      continue
    dp = int(dp_match.group(1))
    mq = int(mq_match.group(1))
    fq = int(fq_match.group(1))
    if dp >= minDepth and dp <= maxDepth and mq >= minMapQ and abs(fq) >= minConsQ:
      mask.addCalledPosition(pos)
      if altAllele != "." and re.match("^[01][/|][01]", genotypes):
        print line,

