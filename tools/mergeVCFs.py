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
      self.file = open(filename, "r")
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
  
    alleles = [ref]
    genotypes = []
    phased = []
    
    for i, l in enumerate(self.current_lines):
      if i not in minIndices:
        genotypes.append([0, 0])
        phased.append(True)
      else:
        al = self.current_lines[i][2]
        genotype = self.current_lines[i][3]
        phased_ = self.current_lines[i][4]
        
        assert al[0] == ref, "inconsistent reference alleles at {}:{}".format(chrom, pos)
        
        new_genotype = []
        for g in genotype:
          if al[g] not in alleles:
            alleles.append(al[g])
          j = alleles.index(al[g])
          new_genotype.append(j)
        genotypes.append(new_genotype)
        phased.append(phased_)
        try:
          self.current_lines[i] = next(self.vcfIterators[i])
        except StopIteration:
          self.current_lines[i] = None
    return (chrom, pos, alleles, genotypes, phased)
  
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
parser.add_argument("files", nargs="+", help="Input VCF and mask files, always given as vcf/mask ordered pair")
parser.add_argument("--sample_names")
args = parser.parse_args()

nrIndidividuals = len(args.files) // 2

if not args.sample_names:
  args.sample_names = ",".join([str(i) for i in range(nrIndidividuals)])

sys.stderr.write("generating msmc input file with {} individuals \n".format(nrIndidividuals))

vcf_files = [f for i, f in enumerate(args.files) if i % 2 == 0]
mask_files = [f for i, f in enumerate(args.files) if i % 2 == 1]

joinedVcfIterator = JoinedVcfIterator(vcf_files)
maskIterators = [MaskIterator(f) for f in mask_files]


# for line in joinedVcfIterator:
#   print(line)
#   print(joinedVcfIterator.current_lines)
#   print(joinedVcfIterator.getMinIndices())

print("##fileformat=VCFv4.1")
print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">')
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format("\t".join(args.sample_names.split(','))))

pos = 0
for chrom, pos, alleles, genotypes, phased in joinedVcfIterator:
  # print(chrom, pos, alleles, genotypes, phased, file=sys.stderr)
  not_missing = False
  for i, m in enumerate(maskIterators):
    if not m.getVal(pos):
      genotypes[i] = (".", ".")
    else:
      not_missing = True
  if not_missing and len(alleles) > 1:
    print("{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT".format(chrom, pos, alleles[0], ",".join(alleles[1:])), end="")
    for g, p in zip(genotypes, phased):
      sep = "|" if p else "/"
      print("\t{}{}{}".format(g[0], sep, g[1]), end="")
    print("")
  
  
  