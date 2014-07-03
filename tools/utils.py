import gzip
import sys
import io

class MaskGenerator:
    def __init__(self, filename, chr):
        self.lastCalledPos = -1
        self.lastStartPos = -1
        self.file = io.TextIOWrapper(gzip.open(filename, "w"))
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

class LegendParser:
    def __init__(self, filename):
        print("opening legend file", filename, file=sys.stderr)
        self.file = io.TextIOWrapper(gzip.open(filename, "r")) if filename[-3:] == ".gz" else open(filename, "r")
        self.end = False
        self.pos = -1
        self.ref_a = ''
        self.alt_a = ''
        self.file.__next__()
    
    def tick(self):
        while True:
            line = self.file.readline().strip()
            if len(line) == 0:
                self.end = True
                break
            fields = line.strip().split()
            if fields[4] == "SNP":
                self.pos = int(fields[1])
                assert fields[2] in "ACTGN" and fields[3] in "ACTG", line
                self.ref_a = fields[2]
                self.alt_a = fields[3]
                break
        
