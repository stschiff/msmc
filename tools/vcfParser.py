#!/usr/bin/env python

import sys
import re
import gzip

class MaskGenerator:
  def __init__(self, filename):
    self.lastCalledPos = -1
    self.lastStartPos = -1
    self.file = gzip.open(filename, "w")
  
  def addCalledPosition(self, pos):
    if self.lastCalledPos == -1:
      self.lastCalledPos = pos
      self.lastStartPos = pos
    elif pos == self.lastCalledPos + 1:
      self.lastCalledPos = pos
    else:
      self.file.write("{}\t{}\n".format(self.lastStartPos, self.lastCalledPos))
      self.lastStartPos = pos
      self.lastCalledPos = pos


class ReferenceSequence:
  def __init__(self, filename, chr):
    self.file = open(filename, "r")
    sys.stderr.write("scanning reference file\n");
    while True:
      self.current_line = self.file.readline().strip()
      if len(self.current_line) == 0:
        sys.stderr.write("reference chromosome not found\n")
        sys.exit(1)
      if not re.match(">" + chr, self.current_line):
        continue
      else:
        sys.stderr.write("found reference file entry\n")
        self.current_line = self.file.readline().strip()
        self.current_pos = 1
        break

  def getBaseAtPos(self, pos):
    if pos < self.current_pos:
      sys.stderr.write("positions not increasing")
      sys.exit(1)
    while pos >= self.current_pos + len(self.current_line):
      self.current_pos += len(self.current_line)
      self.current_line = self.file.readline().strip()
      if len(self.current_line) == 0:
        sys.stderr.write("reference chromosome end reached, position too high\n")
        sys.exit(1)
    return self.current_line[pos - self.current_pos]

minMapQ = 20;
minConsQ = 20;

mean_depth = float(sys.argv[1])
mask_filename = sys.argv[2]
mode = sys.argv[3]
ref = sys.argv[4]
chr = sys.argv[5]

minDepth = mean_depth / 2.0
maxDepth = mean_depth * 2.0

mask = MaskGenerator(mask_filename)

refSeq = None
if mode == "CG":
  refSeq = ReferenceSequence(ref, chr)

lastPos = 0
line_cnt = 0
for line in sys.stdin:
  if line[0] == '#':
    print line,
    continue
  fields = line.strip().split('\t')
  pos = int(fields[1])
  refAllele = fields[3]
  altAllele = fields[4]
  info = fields[7]
  genotypes = fields[9]
  if line_cnt % 10000 == 0:
    sys.stderr.write("parsing position {}\n".format(pos))
  line_cnt += 1

  if mode == "BAM":
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

  elif mode == "CG":
    if pos > lastPos + 1:
      for i in range(lastPos + 1, pos):
        if refSeq.getBaseAtPos(i) != 'N':
          mask.addCalledPosition(i)
    if altAllele == "<CGA_NOCALL>":
      lastPos = int(re.search("END=(\d+)", info).group(1))
    elif len(refAllele) > 1:
      lastPos = pos + len(refAllele) - 1
    elif re.match("^[ACTG]$", altAllele) and re.match("^[01][/|][01]", genotypes) and re.search(":PASS:", genotypes) and refSeq.getBaseAtPos(pos) != 'N':
      print line,
      mask.addCalledPosition(pos)
    lastPos = pos

  else:
    sys.write("unknown mode\n");
    sys.exit(1)
  

