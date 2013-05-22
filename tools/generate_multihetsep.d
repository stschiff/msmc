#!/usr/bin/env rdmd

import std.stdio;
import std.getopt;
import std.functional;
import std.exception;
import std.string;
import std.regex: match;
import std.algorithm;
import std.conv;
import std.array;
import std.math;
import std.range;

char[char[2]] IUPAC_dna;
char[2][char] IUPAC_dna_rev;

enum FILE_TYPE {VCF, CG}

size_t minDepth = 5;
size_t maxDepth = 20;
size_t minMapQ = 20;
size_t minConsQ = 20;
size_t[2][size_t] offspring_pattern;
string[] filenames;
string ref_filename;
FILE_TYPE[] types;
char[][] consensusSequences;
string chromosome;
char[] ref_seq;
string ref_chr;

static this() {
  IUPAC_dna = [
    ['A', 'G'] : 'R',
    ['G', 'A'] : 'R',
    ['C', 'T'] : 'Y',
    ['T', 'C'] : 'Y',
    ['C', 'A'] : 'M',
    ['A', 'C'] : 'M',
    ['T', 'G'] : 'K',
    ['G', 'T'] : 'K',
    ['T', 'A'] : 'W',
    ['A', 'T'] : 'W',
    ['C', 'G'] : 'S',
    ['G', 'C'] : 'S',

    ['A', 'A'] : 'A',
    ['C', 'C'] : 'C',
    ['G', 'G'] : 'G',
    ['T', 'T'] : 'T'
  ];

  IUPAC_dna_rev = [
    'R' : ['A', 'G'],
    'Y' : ['C', 'T'],
    'M' : ['C', 'A'],
    'K' : ['T', 'G'],
    'W' : ['T', 'A'],
    'S' : ['C', 'G'],
    'A' : ['A', 'A'],
    'C' : ['C', 'C'],
    'G' : ['G', 'G'],
    'T' : ['T', 'T']
  ];

}

void main(string[] args) {
  try {
    readArgs(args);
  }
  catch(Exception e) {
    stderr.writeln(e.msg);
    printHelp();
    return;
  }
  
  run();
}

void readArgs(string[] args) {
  getopt(args, std.getopt.config.caseSensitive, std.getopt.config.passThrough, "types|t", toDelegate(&handleTypes), "minDepth|m", &minDepth, "maxDepth|M", &maxDepth, "minMapQ", &minMapQ, "minConsQ", &minConsQ, "offspring|o", toDelegate(&handleOffspringPattern), "chromosome|c", &chromosome, "ref_chr", &ref_chr, "reference_file|r", &ref_filename);
  
  enforce(minDepth > 0 && maxDepth > minDepth, "maxDepth should be larger than minDepth");
  enforce(args.length > 1, "need at least one file");
  filenames = args[1..$];
  if(filenames[0] == "-")
    enforce(filenames.length == 1, "stdin can't be combined with other files");
  enforce(types.length == filenames.length, "types must be set for each file");
  stderr.writeln("input offspring pattern: ", offspring_pattern);
}

void handleTypes(string option, string value) {
  types = value.split(",").map!(a => to!FILE_TYPE(a)).array();
}

void handleOffspringPattern(string option, string value) {
  foreach(pattern; value.split(";")) {
    auto m = pattern.match(r"(\d+)\((\d+),(\d+)\)");
    enforce(m.captures.length() == 4, text(m.captures));
    auto child = m.captures[1].to!size_t();
    auto father = m.captures[2].to!size_t();
    auto mother = m.captures[3].to!size_t();
    offspring_pattern[child] = [father, mother];
  }
}

void printHelp() {
  stderr.writeln("Usage: ./generate_multihetsep.d [options] <file1> [<file2> ...]
if <file1>='-', read from stdin.

General Options:
-t, --types: comma-separated string of either VCF or CG, must be given for each file
-c, --chromosome: chromosome label
-o, --offspring: denote the index of any offspring samples including the two parents, e.g. 1(2,3);4(5,6) to denote two trios where the first and the fourth samples are the offspring of the parental haplotypes (2,3) and (5,6).

Options only for type VCF:
-m, --minDepth [=5]
-M, --maxDepth [=20]
--minMapQ: minimum mapping quality [=20]
--minConsQ: minimum absolute consensus quality [=20]

Options only for type CG:
-r, --reference_file=<file>
--ref_chr=<str>");
}

void run() {
  foreach(fileIndex; 0 .. filenames.length) {
    if(types[fileIndex] == FILE_TYPE.VCF)
      consensusSequences ~= readVCF(filenames[fileIndex]);
    else {
      if(ref_seq.length == 0) {
        readReferenceSequence(ref_filename);
      }
      consensusSequences ~= readCG(filenames[fileIndex]);
    }
  }
  auto maxLength = consensusSequences.map!"a.length"().minCount!"a>b"()[0];
  foreach(c; consensusSequences) {
    c.length = maxLength;
  }
  printOutput();
}

char[][] readVCF(string filename) {
  char[][] ret;
  auto nrSamples = 0UL;
  auto file = stdin;
  if(filename != "-") {
    stderr.writeln("reading file ", filename);
    file = File(filename, "r");
  }
  else {
    stderr.writeln("reading from stdin");
  }
  auto cnt = 0;
  foreach(line; file.byLine().filter!(l => l.startsWith(chromosome))()) {
    auto fields = line.strip().split("\t");
    if(nrSamples == 0) {
      nrSamples = fields.length - 9;
      ret = new char[][nrSamples];
    }
    auto pos = fields[1].to!size_t();
    if(cnt++ % 100000 == 0)
      stderr.writeln("scanning pos ", pos);
    auto i = pos - 1;
    if(i >= ret[0].length) {
      auto previousLength = ret[0].length;
      foreach(j; 0 .. nrSamples) {
        ret[j].length = i + 1;
        ret[j][previousLength .. $] = 'N';
      }
    }
    if(fields[3].match(r"^[ACTGactg]$") && fields[4].match(r"^[ACTGactg\.]$")) {
      auto dp = fields[7].match(r"DP=(\d+)").captures[1].to!size_t();
      auto mq = fields[7].match(r"MQ=(\d+)").captures[1].to!int();
      auto fq = fields[7].match(r"FQ=([\d-]+)").captures[1].to!int();
      if(dp >= minDepth && dp <= maxDepth && mq >= minMapQ && abs(fq) >= minConsQ) {
        if(fields[4] != ".") {
          foreach(j; 0 .. nrSamples) {
            auto gen1 = [fields[9 + j][0]].to!int();
            auto gen2 = [fields[9 + j][2]].to!int();
            char[2] dinuc = [fields[3 + gen1][0], fields[3 + gen2][0]];
            assert(dinuc in IUPAC_dna);
            ret[j][i] = IUPAC_dna[dinuc];
          }
        }
        else {
          foreach(j; 0 .. nrSamples) {
            ret[j][i] = fields[3][0];
          }
        }
      }
    }
  }
  return ret;
}

char[] readCG(string filename) {
  enum zygosity_t {NO_CALL, HAP, HALF, HOM, HET_REF, HET_ALT}
  enum vartype_t {SNP, INS, DEL, SUB, REF, COMPLEX, NO_REF, PAR_CALLED_IN_X}
  enum varquality_t {VQLOW, VQHIGH}
  
  stderr.writeln("reading file ", filename);
  auto cg_file = File(filename, "r");
  auto line_count = 0;
  auto seq = new char[ref_seq.length];
  seq[] = 'N';
  foreach(line; cg_file.byLine().filter!(l => !(startsWith(l, "#") || startsWith(l, ">") || l.strip().length == 0))()) {
    if(line_count++ % 10000 == 0)
      stderr.writeln("processing line ", line_count);

    auto fields = line.strip().split("\t");
    if(fields[2] != chromosome) continue;

    auto chrom = fields[2];
    auto begin = to!int(fields[3]);
    auto end = to!int(fields[4]);
    auto zygosity = to!zygosity_t(tr(fields[5], "-", "_").toUpper());
    auto varType = to!vartype_t(tr(fields[6], "-", "_").toUpper());
    auto alleleRef = fields[7];
    auto allele1 = fields[8];
    auto allele2 = fields[9];
    
    if(varType == vartype_t.REF && zygosity == zygosity_t.HOM) {
      seq[begin .. end] = ref_seq[begin .. end];
    }
    if(varType == vartype_t.SNP) {
      if((zygosity == zygosity_t.HOM || zygosity == zygosity_t.HET_REF || 
         zygosity == zygosity_t.HET_ALT) && ref_seq[begin] != 'N')
      {
        assert(end - begin == 1, text(fields));
        assert(
          alleleRef[0] == ref_seq[begin],
          format("%s %s cg_ref:%s my_ref:%s", chrom, end, alleleRef[0], ref_seq[begin])
        );
        auto allele1_qual = to!varquality_t(fields[14]);
        auto allele2_qual = to!varquality_t(fields[15]);
        if(allele1_qual == varquality_t.VQHIGH &&
           allele2_qual == varquality_t.VQHIGH)
        {
          char[2] dinuc = [allele1[0], allele2[0]];
          assert(dinuc in IUPAC_dna);
          seq[begin] = IUPAC_dna[dinuc];
        }
        // stderr.writefln("%s cons:%s ref:%s", end, cg_sequences[cg_file_index][begin], ref_seq[begin]);
      }
    }
  }
  return seq;
}

void readReferenceSequence(string filename) {
  auto ref_file = File(filename, "r");
  auto tag = ">" ~ ref_chr ~ " ";
  char[] line;
  do
    ref_file.readln(line);
  while(!startsWith(line, tag) && !ref_file.eof());
  
  stderr.writeln("found fasta region ", strip(line));
  ref_file.readln(line);
  do {
    ref_seq ~= line.strip();
    ref_file.readln(line);
  } while(!startsWith(line, ">") && !ref_file.eof());

  auto L = ref_seq.length;
  stderr.writeln("read ", L, " nucleotides");
}

void printOutput() {
  auto nr_called_sites = 0;
  auto nrHaplotypes = consensusSequences.length;
  foreach(i; 0 .. consensusSequences[0].length) {
    if(i % 100000 == 0)
      stderr.writeln("position ", i);
    if(!hasGap(i)) {
      nr_called_sites += 1;
      if(hasSNP(i)) {
        auto allPhasings = getAllPhasings(consensusSequences.map!(c => c[i])().array());
        auto nrHaplotypes = allPhasings[0].length;
        auto prunedPhasings = pruneInconsistentPhasings(allPhasings, offspring_pattern);
        if(prunedPhasings.length == 0)
          prunedPhasings = [iota(nrHaplotypes).map!"'?'"().array()];
        auto alleleString = prunedPhasings.joiner(",").array();
        writefln("%s %s %s %s", chromosome, i + 1, nr_called_sites, alleleString);
        nr_called_sites = 0;
      }
    }
  }
}

bool hasGap(size_t i) {
  foreach(s; consensusSequences) {
    if(s[i] !in IUPAC_dna_rev)
      return true;
  }
  return false;
}

bool hasSNP(size_t i) {
  foreach(s; consensusSequences) {
    if(!canFind("ACTG", s[i]) || s[i] != consensusSequences[0][i])
      return true;
  }
  return false;
}

char[][] getAllPhasings(in char[] genotypes) {
  char[][] ret;
  auto alleles = IUPAC_dna_rev[genotypes[0]];
  if(alleles[0] == alleles[1])
    ret = [alleles.dup];
  else {
    ret = [alleles.dup, [alleles[1], alleles[0]]];
  }
  if(genotypes.length > 1) {
    foreach(g; genotypes[1 .. $]) {
      auto al = IUPAC_dna_rev[g];
      if(al[0] == al[1]) {
        foreach(ref str; ret) {
          str ~= al.dup;
        }
      }
      else {
        char[][] newRet;
        foreach(r; ret) {
          newRet ~= (r ~ al.dup);
          newRet ~= (r ~ [al[1], al[0]]);
        }
        ret = newRet;
      }
    }
  }
  return ret.uniq.array();
}

unittest {
  assert(getAllPhasings(['A']) == ["AA"]);
  assert(getAllPhasings(['A', 'C']) == ["AACC"], text(getAllPhasings(['A', 'C'])));
  assert(getAllPhasings(['R']) == ["AG", "GA"]);
  assert(getAllPhasings(['R', 'A']) == ["AGAA", "GAAA"]);
  assert(getAllPhasings(['R', 'R']) == ["AGAG", "AGGA", "GAAG", "GAGA"]);
}

char[][] pruneInconsistentPhasings(in char[][] allPhasings, size_t[2][size_t] offspring_pattern) {
  size_t[] inconsistentIndices;
  foreach(offspringIndex; offspring_pattern.keys()) {
    foreach(i, alleles; allPhasings) {
      if(!consistentPhasing(alleles, offspring_pattern[offspringIndex])) {
        inconsistentIndices ~= i;
      }
    }
  }
  
  auto consistentPhasings = iota(allPhasings.length)
         .filter!(i => !canFind(inconsistentIndices, i))()
         .map!(i => allPhasings[i])();
  return consistentPhasings.map!(al => iota(al.length).filter!(i => i !in offspring_pattern)().array()).array();
}

unittest {
  
}
    