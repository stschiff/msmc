#!/usr/bin/env rdmd
  
import std.stdio;
import std.getopt;
import std.exception;
import std.file;
import std.conv;
import std.string;
import std.algorithm;
import std.array;
import std.range;

string haplotypeFileName;
string multihetsep_filename;
size_t[] multihetsepIndices;
size_t[] hapsIndices;
size_t nrIndividuals;

void main(string[] args) {
  try {
    readArgs(args);
  }
  catch(Exception e) {
    stderr.writeln(e.msg);
    printHelpMessage();
    return;
  }
  run();
} 

void readArgs(string[] args) {
  getopt(args, "hapsFile|h", &haplotypeFileName, "multihetsepIndex|i", &multihetsepIndices, "hapsIndex|j", &hapsIndices);
  nrIndividuals = multihetsepIndices.length;
  enforce(nrIndividuals > 0, "need to specify at least one individual");
  enforce(multihetsepIndices.length == hapsIndices.length, "need same number of indices in multihetsep and hapsfile");
  enforce(isFile(haplotypeFileName), text("file ", haplotypeFileName, " not found"));
  enforce(args.length == 2, "need one multihetsep file");
  multihetsep_filename = args[1];
  enforce(isFile(multihetsep_filename), text("file ", multihetsep_filename, " not found"));
}


void printHelpMessage() {
  stderr.writeln("./multihetsep_phase.d [Options] <multihetsep_file>
Options:
  -h, --hapsFile
  -i, --multihetsepIndex
  -j, --hapsIndex");
}

void run() {
  auto lastHapsFilePos = 0UL;
  auto multihetsepF = File(multihetsep_filename, "r");
  auto hapsFile = File(haplotypeFileName, "r");
  auto hapsFileRange = hapsFile.byLine;
  char[][] hapsFields;
  foreach(line; multihetsepF.byLine) {
    // stderr.writeln(line);
    auto fields = line.strip.split;
    auto chr = fields[0];
    auto pos = fields[1].to!size_t();
    auto nr_called = fields[2];
    while(lastHapsFilePos < pos && !hapsFileRange.empty) {
      auto hapsFileLine = hapsFileRange.front.dup;
      hapsFileRange.popFront;
      hapsFields = hapsFileLine.strip.split;
      lastHapsFilePos = hapsFields[2].to!size_t();
    }
    // stderr.writeln(hapsFields);
    char[] alleleString;
    if(lastHapsFilePos == pos) {
      foreach(i; 0 .. nrIndividuals) {
        auto alleles = fields[3 + multihetsepIndices[i]];
        auto al = hapsFields[3..5].joiner.array;
        // stderr.writeln(hapsFields[5 + 2 * hapsIndices[i] .. 5 + 2 * (hapsIndices[i] + 1)]);
        auto gens = hapsFields[5 + 2 * hapsIndices[i] .. 5 + 2 * (hapsIndices[i] + 1)].map!(a => a.to!size_t())().array;
        enforce((alleles[0] == al[0] && alleles[1] == al[1]) || (alleles[0] == al[1] || alleles[1] == al[0]), 
                text("wrong alleles at position", pos));
        // stderr.writeln(pos, " ", gens, " ", al);
        alleleString ~= gens.map!(g => al[g].to!char())().array;
      }
    }
    else {
      char[][] phasings;
      foreach(i; 0 .. nrIndividuals) {
        auto alleles = fields[3 + multihetsepIndices[i]];
        if(phasings.length == 0) {
          if(alleles[0] == alleles[1])
            phasings = [alleles.dup];
          else {
            phasings = [alleles.dup, alleles.dup.retro.map!"a.to!char()"().array];
          }
        }
        else {
          if(alleles[0] == alleles[1]) {
            foreach(ref str; phasings) {
              str ~= alleles.dup;
            }
          }
          else {
            char[][] newPhasings;
            foreach(r; phasings) {
              newPhasings ~= (r ~ alleles.dup);
              newPhasings ~= (r ~ alleles.dup.retro.map!"a.to!char()"().array);
            }
            phasings = newPhasings;
          }
        }
      }
      alleleString = phasings.joiner(",").map!"a.to!char()"().array;
    }
    writefln("%s\t%s\t%s\t%s", chr, pos, nr_called, alleleString);
  }
}
