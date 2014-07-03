#!/usr/bin/env python3

import sys
import gzip
import string
import copy
import argparse
import io

class MaskIterator:
  def __init__(self, filename, negative=False):
    if filename[-3:] == ".gz":
      self.file = io.TextIOWrapper(gzip.open(filename, "r"))
    else:
      self.file = io.TextIOWrapper(open(filename, "r"))
    self.eof = False
    self.lastPos = 1
    self.negative = negative
    self.readLine()

  def readLine(self):
    try:
      line = next(self.file)
      fields = line.strip().split()
      if len(fields) == 2:
        self.start = int(fields[0])
        self.end = int(fields[1])
      else:
        self.start = int(fields[1]) + 1
        self.end = int(fields[2])
    except StopIteration:
      self.eof = True
  
  def getVal(self, pos):
    assert pos >= self.lastPos
    self.lastPos = pos
    while pos > self.end and not self.eof:
      self.readLine()
    if pos >= self.start and pos <= self.end:
      return True if not self.negative else False
    else:
      return False if not self.negative else True

class MergedMask:
  def __init__(self, mask_iterators):
    self.maskIterators = mask_iterators

  def getVal(self, pos):
    return all((m.getVal(pos) for m in self.maskIterators))

class VcfIterator:
  def __init__(self, filename):
    self.file = io.TextIOWrapper(gzip.open(filename, "r"))
  
  def __iter__(self):
    return self
  
  def __next__(self):
    line = next(self.file)
    while line[0] == "#":
      line = next(self.file)
    fields = line.strip().split()
    chrom = fields[0]
    pos = int(fields[1])
    alleles = [fields[3]]
    for alt_a in fields[4].split(","):
      alleles.append(alt_a)
    geno = fields[9][:3]
    phased = geno[1] == "|"
    return (chrom, pos, tuple(alleles), (int(geno[0]), int(geno[2])), phased)

class OrderedAlleles:
  def __init__(self):
    self.ordered_alleles = []
  
  def addGenotype(self, a1, a2, phasing):
    if len(self.ordered_alleles) == 0:
      self.ordered_alleles = [[a1, a2]]
      if not phasing and a1 != a2:
        self.ordered_alleles.append([a2, a1])
    else:
      new = []
      for o in self.ordered_alleles:
        new.append(o + [a1, a2])
        if not phasing and a1 != a2:
          new.append(o + [a2, a1])
      self.ordered_alleles = new

  def getPrint(self):
    if len(self.ordered_alleles[0]) == 2 or len(self.ordered_alleles) == 1:
      return ''.join(self.ordered_alleles[0])
    else:
      return ','.join([''.join(o) for o in self.ordered_alleles])

class JoinedVcfIterator:
  def __init__(self, filenames):
    self.vcfIterators = [VcfIterator(f) for f in filenames]
    self.current_lines = [next(v) for v in self.vcfIterators]
  
  def __iter__(self):
    return self
  
  def __next__(self):
    minIndices = self.getMinIndices()
    chrom = self.current_lines[minIndices[0]][0]
    pos = self.current_lines[minIndices[0]][1]
    ref = self.current_lines[minIndices[0]][2][0]
  
    ordered_alleles = OrderedAlleles()
    
    for i, l in enumerate(self.current_lines):
      if i not in minIndices:
        ordered_alleles.addGenotype(ref, ref, True)
      else:
        alleles = self.current_lines[i][2]
        geno = self.current_lines[i][3]
        phased = self.current_lines[i][4]
        ordered_alleles.addGenotype(alleles[geno[0]], alleles[geno[1]], phased)
        try:
          self.current_lines[i] = next(self.vcfIterators[i])
        except StopIteration:
          self.current_lines[i] = None
    return (chrom, pos, ordered_alleles.getPrint())
  
  def getMinIndices(self):
    activeLines = [(i, l) for i, l in enumerate(self.current_lines) if l]
    if len(activeLines) == 0:
      raise StopIteration
    if len(activeLines) == 1:
      return [activeLines[0][0]]
    else:
      minIndices = [activeLines[0][0]]
      minPos = activeLines[0][1][1]
      for a in activeLines[1:]:
        if a[1][1] == minPos:
          minIndices.append(a[0])
        if a[1][1] < minPos:
          minPos = a[1][1]
          minIndices = [a[0]]
      return minIndices
    

parser = argparse.ArgumentParser()
parser.add_argument("files", nargs="+", help="Input VCF files")
parser.add_argument("--mask", action="append", help="apply masks in bed format, should be given once for the calling mask from each individual, and in addition can be given for e.g. mappability or admixture masks")
parser.add_argument("--negative_mask", action="append", help="same as mask, but interpreted as negative mask, so places where sites should be excluded")
args = parser.parse_args()

nrIndidividuals = len(args.files)
nrHaplotypes = 2 * nrIndidividuals

sys.stderr.write("generating msmc input file with {} haplotypes\n".format(nrHaplotypes))

joinedVcfIterator = JoinedVcfIterator(args.files)
maskIterators = []
if args.mask:
  for f in args.mask:
    sys.stderr.write("adding mask: {}\n".format(f))
    maskIterators.append(MaskIterator(f))
if args.negative_mask:
  for nm in args.negative_mask:
    sys.stderr.write("adding negative mask: {}\n".format(nm))
    maskIterators.append(MaskIterator(nm, True))

mergedMask = MergedMask(maskIterators)

def is_segregating(alleles):
  orders = alleles.split(",")
  for o in orders:
    for a in o[1:]:
      if a != o[0]:
        return True
  return False

pos = 0
nr_called = 0
for chrom, snp_pos, alleles in joinedVcfIterator:
  # sys.stderr.write("{}\t{}\t{}\n".format(chrom, snp_pos, alleles))
  while pos < snp_pos:
    pos += 1
    if mergedMask.getVal(pos):
      nr_called += 1
    if pos % 1000000 == 0:
      print("processing pos {}".format(pos), file=sys.stderr)
  if mergedMask.getVal(snp_pos):
    if is_segregating(alleles):
      print(chrom, snp_pos, nr_called, alleles, sep="\t")
      nr_called = 0
  
  
  