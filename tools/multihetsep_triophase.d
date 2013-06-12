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

size_t[2][size_t] offspring_pattern;
string filename;

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
  getopt(args, std.getopt.config.caseSensitive, std.getopt.config.passThrough, "offspring|o", toDelegate(&handleOffspringPattern));
  
  enforce(args.length == 2, "need one file");
  filename = args[1];
  stderr.writeln("input offspring pattern: ", offspring_pattern);
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
  stderr.writeln("Usage: ./multihetsep_triophase.d [options] <multihetsep_file>

Options:
-o, --offspring: denote the index of any offspring samples including the two parents, e.g. 1(2,3);4(5,6) to denote two trios where the first and the fourth samples are the offspring of the parental haplotypes (2,3) and (5,6).");
}

void run() {
  auto f = File(filename, "r");
  
  foreach(line; f.byLine) {
    auto fields = line.strip.split;
    auto genotypes = fields[3 .. $];
    auto nrHaplotypes = 2 * genotypes.length;
    auto allPhasings = getAllPhasings(genotypes);
    char[] alleleString;
    if(nrHaplotypes == 2)
      alleleString = allPhasings[0];
    else {
      auto prunedPhasings = pruneInconsistentPhasings(allPhasings, offspring_pattern);
      if(prunedPhasings.length == 0)
        prunedPhasings = [iota(nrHaplotypes - offspring_pattern.length * 2).map!"'?'"().array()];
      alleleString = prunedPhasings.joiner(",").array().to!(char[]);
    }
    writefln("%s\t%s\t%s\t%s", fields[0], fields[1], fields[2], alleleString);
  }
}

char[][] getAllPhasings(in char[][] genotypes) {
  char[][] ret;
  if(genotypes[0][0] == genotypes[0][1])
    ret = [genotypes[0].dup];
  else {
    ret = [genotypes[0].dup, genotypes[0].dup.retro.map!"a.to!char()"().array];
  }
  if(genotypes.length > 1) {
    foreach(g; genotypes[1 .. $]) {
      if(g[0] == g[1]) {
        foreach(ref str; ret) {
          str ~= g.dup;
        }
      }
      else {
        char[][] newRet;
        foreach(r; ret) {
          newRet ~= (r ~ g.dup);
          newRet ~= (r ~ g.dup.retro.map!"a.to!char()"().array);
        }
        ret = newRet;
      }
    }
  }
  return ret.uniq.array();
}

unittest {
  assert(equal(getAllPhasings([['A', 'A']]).joiner(","), "AA"));
  assert(equal(getAllPhasings([['A', 'A'], ['C', 'C']]).joiner(","), "AACC"));
  assert(equal(getAllPhasings([['A', 'G']]).joiner(","), "AG,GA"));
  assert(equal(getAllPhasings([['A', 'G'], ['A', 'A']]).joiner(","), "AGAA,GAAA"));
  assert(equal(getAllPhasings([['A', 'G'], ['A', 'G']]).joiner(","), "AGAG,AGGA,GAAG,GAGA"));
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