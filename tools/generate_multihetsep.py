#!/usr/bin/env python

import sys
import gzip
import string
import copy

if len(sys.argv) < 2:
  sys.stderr.write("""Usage: ./generate_multihetsep.py <vcf_1> <mask_1> <vcf_2> <mask_2> ...
where both vcf files and mask files must be gzip-compressed
""")

filenames = sys.argv[1:]
nrIndidividuals = len(filenames) / 2
nrHaplotypes = 2 * nrIndidividuals

class MaskIterator:
  def __init__(self, filename):
    self.file = gzip.open(filename, "r")
    self.eof = False
    self.lastPos = 1
    self.readLine()

  def readLine(self):
    try:
      line = self.file.next()
      fields = string.split(string.strip(line))
      self.start = int(fields[0])
      self.end = int(fields[1])
    except StopIteration:
      self.eof = True
  
  def getVal(self, pos):
    assert pos >= self.lastPos
    self.lastPos = pos
    while pos > self.end and not self.eof:
      self.readLine()
    if pos >= self.start and pos <= self.end:
      return True
    else:
      return False

class MergedMask:
  def __init__(self, mask_files):
    self.maskIterators = [MaskIterator(f) for f in mask_files]

  def getVal(self, pos):
    return all((m.getVal(pos) for m in self.maskIterators))

class VcfIterator:
  def __init__(self, filename):
    self.file = gzip.open(filename, "r")
  
  def __iter__(self):
    return self
  
  def next(self):
    line = self.file.next()
    while line[0] == "#":
      line = self.file.next()
    fields = string.split(string.strip(line))
    chrom = fields[0]
    pos = int(fields[1])
    alleles = fields[3:5]
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
    self.current_lines = [v.next() for v in self.vcfIterators]
  
  def __iter__(self):
    return self
  
  def next(self):
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
          self.current_lines[i] = self.vcfIterators[i].next()
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
    
    
# DIR="/lustre/scratch113/projects/msmc/Ethiopians/consensus_calls"
# joined = JoinedVcfIterator([DIR + "/NA12878/NA12878_chr1.phased.vcf.gz", DIR + "/NA19240/NA19240_chr1.phased.vcf.gz", DIR + "/egpg5305762/egpg5305762_chr1.phased.vcf.gz", DIR + "/egpg5305814/egpg5305814_chr1.phased.vcf.gz"])
# 
# for line in joined:
#   print line


sys.stderr.write("generating msmc input file with {} haplotypes\n".format(nrHaplotypes))

joinedVcfIterator = JoinedVcfIterator([filenames[2 * i] for i in range(nrIndidividuals)])
mergedMask = MergedMask([filenames[2 * i + 1] for i in range(nrIndidividuals)])

pos = 0
nr_called = 0
for chrom, snp_pos, alleles in joinedVcfIterator:
  while pos < snp_pos:
    pos += 1
    if mergedMask.getVal(pos):
      nr_called += 1
    if pos % 1000000 == 0:
      sys.stderr.write("processing pos {}\n".format(pos))
  if mergedMask.getVal(snp_pos):
    print "{}\t{}\t{}\t{}".format(chrom, snp_pos, nr_called, alleles)
    nr_called = 0
  
  
  