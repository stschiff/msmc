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
string chromosome;
string ref_chr;
bool diploidOutput;
string positions_filename;

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
    ['T', 'T'] : 'T',
      
    ['N', 'N'] : 'N'
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
    'T' : ['T', 'T'],
    'N' : ['N', 'N']
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
  getopt(args, std.getopt.config.caseSensitive, std.getopt.config.passThrough, "types|t", toDelegate(&handleTypes), "minDepth|m", &minDepth, "maxDepth|M", &maxDepth, "minMapQ", &minMapQ, "minConsQ", &minConsQ, "offspring|o", toDelegate(&handleOffspringPattern), "chromosome|c", &chromosome, "ref_chr", &ref_chr, "reference_file|r", &ref_filename, "diploidOutput|d", &diploidOutput, "positions_file|p", &positions_filename);
  
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
-d, --diploidOutput: Output as diploid alleles, using IUPAC symbols for heterozygotes
-p, --positions_file: File with line-separated positions, at which output is forced (even for missing or homozygous calls)

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
  char[] ref_seq;
  char[][] consensusSequences;
  foreach(fileIndex; 0 .. filenames.length) {
    if(types[fileIndex] == FILE_TYPE.VCF)
      consensusSequences ~= readVCF(filenames[fileIndex]);
    else {
      if(ref_seq.length == 0) {
        ref_seq = readFastaSequence(ref_filename, ref_chr);
      }
      consensusSequences ~= readCG(filenames[fileIndex], ref_seq);
    }
  }
  normalizeLengths(consensusSequences, 'N');
  printOutput(consensusSequences);
}

char[][] readVCF(string filename) {
  char[][] ret;
  auto nrSamples = 0UL;
  auto file = openFile(filename);
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
    
    foreach(j; 0 .. nrSamples)
      fillLength(ret[j], i + 1, 'N');
    
    if(fields[3].match(r"^[ACTGactg]$") && fields[4].match(r"^[ACTGactg\.]$")) {
      auto dp_match = fields[7].match(r"DP=(\d+)");
      auto mq_match = fields[7].match(r"MQ=(\d+)");
      auto fq_match = fields[7].match(r"FQ=([\d-]+)");
      if(!(dp_match && mq_match && fq_match)) // this bypasses a bug in samtools that causes a lack of fields in the INFO column.
        continue;
      auto dp = dp_match.captures[1].to!size_t();
      auto mq = mq_match.captures[1].to!int();
      auto fq = fq_match.captures[1].to!int();
            
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

File openFile(string filename) {
  auto file = stdin;
  if(filename != "-") {
    stderr.writeln("reading file ", filename);
    file = File(filename, "r");
  }
  else {
    stderr.writeln("reading from stdin");
  }
  return file;
}

void fillLength(ref char[] seq, size_t newLength, char fill) {
  if(newLength > seq.length) {
    auto previousLength = seq.length;
    seq.length = newLength;
    seq[previousLength .. $] = fill;
  }
}

unittest {
  auto seq = "ACCT".dup;
  fillLength(seq, 8UL, 'N');
  assert(equal(seq, "ACCTNNNN"), text(seq));
  fillLength(seq, 6UL, 'N');
  assert(equal(seq, "ACCTNNNN"));
}

char[] readCG(string filename, char[] ref_seq) {
  enum zygosity_t {NO_CALL, HAP, HALF, HOM, HET_REF, HET_ALT}
  enum vartype_t {SNP, INS, DEL, SUB, REF, COMPLEX, NO_REF, PAR_CALLED_IN_X}
  enum varquality_t {VQLOW, VQHIGH}
  
  auto cg_file = openFile(filename);
  auto line_count = 0;
  auto seq = new char[ref_seq.length];
  seq[] = 'N';
  bool chromosome_read;
  foreach(line; cg_file.byLine().filter!(l => !(startsWith(l, "#") || startsWith(l, ">") || l.strip().length == 0))()) {

    auto fields = line.strip().split("\t");

    auto chrom = fields[2];
    auto begin = to!int(fields[3]);
    auto end = to!int(fields[4]);
    if(line_count++ % 100000 == 0) {
      stderr.writefln("processing pos %s:%s-%s", chrom, begin, end);
      // break;
    }
    if(chrom != chromosome) {
      if(chromosome_read)
        break;
      else
        continue;
    }
    chromosome_read = true;
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

char[] readFastaSequence(string filename, string ref_chr) {
  auto ref_file = File(filename, "r");
  auto tag = ">" ~ ref_chr ~ " ";
  char[] line;
  do
    ref_file.readln(line);
  while(!startsWith(line, tag) && !ref_file.eof());
  
  stderr.writeln("found fasta region ", strip(line));
  ref_file.readln(line);
  char[] ref_seq;
  do {
    ref_seq ~= line.strip();
    ref_file.readln(line);
  } while(!startsWith(line, ">") && !ref_file.eof());

  auto L = ref_seq.length;
  stderr.writeln("read ", L, " nucleotides");
  return ref_seq;
}

void normalizeLengths(char[][] consensusSequences, char fill) {
  auto maxLength = consensusSequences.map!"a.length"().minCount!"a>b"()[0];
  stderr.writeln("normalizing to length ", maxLength);
  foreach(ref c; consensusSequences) {
    fillLength(c, maxLength, fill);
  }
}

unittest {
  auto seqs = ["ACCT", "ACCCTTT", "ACA"].map!"a.dup"().array();
  char fill = 'N';
  normalizeLengths(seqs, fill);
  assert(equal(seqs[0], "ACCTNNN"), text(seqs[0]));
  assert(equal(seqs[1], "ACCCTTT"));
  assert(equal(seqs[2], "ACANNNN"));
}

void printOutput(char[][] consensusSequences) {
  auto nr_called_sites = 0;
  auto nrHaplotypes = 2 * consensusSequences.length;
  size_t[] forced_positions;
  size_t current_position_index = 0;
  if(positions_filename.length > 0) {
    forced_positions = readPositions(positions_filename);
    stderr.writeln("read positions: ", forced_positions[0..20]);
  }
  foreach(i; 0 .. consensusSequences[0].length) {
    if(i % 1000000 == 0)
      stderr.writeln("position ", i);
    if(hasGap(consensusSequences, i))
      continue;
    nr_called_sites += 1;
    if(hasSNP(consensusSequences, i) || ((current_position_index < forced_positions.length - 1) && (forced_positions[current_position_index] == i + 1))) {
      char[] alleleString;
      if(diploidOutput) {
        alleleString = consensusSequences.map!(c => c[i])().array();
      }
      else {
        auto allPhasings = getAllPhasings(consensusSequences.map!(c => c[i])().array());
        if(nrHaplotypes == 2)
          alleleString = allPhasings[0];
        else {
          auto prunedPhasings = pruneInconsistentPhasings(allPhasings, offspring_pattern);
          if(prunedPhasings.length == 0)
            prunedPhasings = [iota(nrHaplotypes - offspring_pattern.length * 2).map!"'?'"().array()];
          alleleString = prunedPhasings.joiner(",").array().to!(char[]);
        }
      }
      if(forced_positions[current_position_index] == i + 1) {
        // writeln("forced position: ", forced_positions[current_position_index]);
        current_position_index += 1;
      }
      writefln("%s\t%s\t%s\t%s", chromosome, i + 1, nr_called_sites, alleleString);
      nr_called_sites = 0;
    }
  }
}

size_t[] readPositions(string filename) {
  auto f = File(filename, "r");
  return f.byLine().map!"to!size_t(a.strip())"().array();
}

bool hasGap(in char[][] consensusSequences, size_t i) {
  foreach(s; consensusSequences) {
    if(s[i] !in IUPAC_dna_rev || s[i] == 'N')
      return true;
  }
  return false;
}

unittest {
  auto seqs = ["ACCC", "ANCC"].map!"a.dup"().array();
  assert(hasGap(seqs, 1));
  assert(!hasGap(seqs, 0));
}

bool hasSNP(in char[][] consensusSequences, size_t i) {
  foreach(s; consensusSequences) {
    if(!canFind("ACTGN", s[i]) || s[i] != consensusSequences[0][i])
      return true;
  }
  return false;
}

unittest {
  auto seqs = ["ACCR", "ANCC"].map!"a.dup"().array();
  assert(!hasSNP(seqs, 0));
  assert(hasSNP(seqs, 3));
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
  assert(equal(getAllPhasings(['A']).joiner(","), "AA"));
  assert(equal(getAllPhasings(['A', 'C']).joiner(","), "AACC"));
  assert(equal(getAllPhasings(['R']).joiner(","), "AG,GA"));
  assert(equal(getAllPhasings(['R', 'A']).joiner(","), "AGAA,GAAA"));
  assert(equal(getAllPhasings(['R', 'R']).joiner(","), "AGAG,AGGA,GAAG,GAGA"));
}

char[][] pruneInconsistentPhasings(in char[][] allPhasings, size_t[2][size_t] offspring_pattern) {
  size_t[] consistentIndices;
  foreach(i, alleles; allPhasings) {
    auto consistent = true;
    foreach(offspringIndex; offspring_pattern.keys()) {
      if(!consistentPhasing(alleles, offspringIndex, offspring_pattern[offspringIndex])) {
        consistent = false;
        break;
      }
    }
    if(consistent)
      consistentIndices ~= i;
  }
  
  return consistentIndices
    .map!(i => allPhasings[i].dup)()
    .map!(al => iota(al.length).filter!(i => (i/2 + 1) !in offspring_pattern)().map!(i => al[i]).array())
    .array();
}

unittest {
  auto seqs = ["ACAAAC", "ACAACA", "CAAAAC", "CAAACA"].map!"a.dup"().array();
  size_t[2][size_t] o;
  o[3] = [1, 2];
  auto pruned = pruneInconsistentPhasings(seqs, o);
  assert(equal(pruned, ["CAAA"]), text(pruned));
}

bool consistentPhasing(in char[] alleles, size_t offspring, size_t[2] parents) {
  auto offspringAl1 = (offspring - 1) * 2;
  auto offspringAl2 = (offspring - 1) * 2 + 1;
  auto fatherIndex = (parents[0] - 1) * 2;
  auto motherIndex = (parents[1] - 1) * 2;
  return alleles[offspringAl1] == alleles[fatherIndex] && alleles[offspringAl2] == alleles[motherIndex];
}

unittest {
  // parents hom, child het
  assert(consistentPhasing("AACCAC", 3, [1, 2]));
  assert(!consistentPhasing("AACCCA", 3, [1, 2]));
  
  // one parent het
  assert(consistentPhasing("ACAAAA", 3, [1, 2]));
  assert(!consistentPhasing("CAAAAA", 3, [1, 2]));
  assert(consistentPhasing("AAACAA", 3, [1, 2]));
  assert(!consistentPhasing("AACAAA", 3, [1, 2]));

  // both parents het
  assert(consistentPhasing("ACACAA", 3, [1, 2]));
  assert(!consistentPhasing("ACCAAA", 3, [1, 2]));
  assert(!consistentPhasing("CAACAA", 3, [1, 2]));
  assert(!consistentPhasing("CACAAA", 3, [1, 2]));
  
  // all hets
  assert(consistentPhasing("ACCAAC", 3, [1, 2]));
  assert(consistentPhasing("CAACCA", 3, [1, 2]));
  assert(!consistentPhasing("ACACAC", 3, [1, 2]));
  assert(!consistentPhasing("CACAAC", 3, [1, 2]));
  assert(!consistentPhasing("ACACCA", 3, [1, 2]));
  assert(!consistentPhasing("CACACA", 3, [1, 2]));
  
}