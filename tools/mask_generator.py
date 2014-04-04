import gzip

class MaskGenerator:
  def __init__(self, filename, chr):
    self.lastCalledPos = -1
    self.lastStartPos = -1
    self.file = gzip.open(filename, "w")
    self.chr = chr
  
  # assume 1-based coordinate, output in bed format
  def addCalledPosition(self, pos):
    if self.lastCalledPos == -1:
      self.lastCalledPos = pos
      self.lastStartPos = pos
    elif pos == self.lastCalledPos + 1:
      self.lastCalledPos = pos
    else:
      self.file.write("{}\t{}\t{}\n".format(self.chr, self.lastStartPos - 1, self.lastCalledPos))
      self.lastStartPos = pos
      self.lastCalledPos = pos


