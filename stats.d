import std.stdio;
import std.getopt;
import std.string;
import std.conv;
import std.array;
import std.algorithm;
import std.c.stdlib;
import std.exception;
import utils;
import data;
import time_intervals;


class StatsApplication {
  
  string[] fileNames;
  size_t nrHaplotypes;
  
  static void runWithArgs(string[] args) {
    auto statsApplication = new StatsApplication(args);
    statsApplication.run();
  }
  
  this(string[] args) {
    parseCommandLine(args);
  }
  
  void parseCommandLine(string[] args) {
    if(args.length == 1)
      displayHelpMessageAndExit();
    try {
      readArguments(args);
      validate();
    }
    catch(Exception e) {
      displayHelpMessageAndExit(e);
    }
  }

  void readArguments(string[] args) {
    fileNames = args[1..$];
    nrHaplotypes = getNrHaplotypesFromFile(fileNames[0]);
    stderr.writeln("found ", nrHaplotypes, " haplotypes in file");
  }

  private void validate() {
    enforce(fileNames.length > 0, "need at least one input file");
    enforce(nrHaplotypes >= 2, "need at least two haplotypes");
  }
  
  
  static void displayHelpMessageAndExit() {
    stderr.writeln("msmc stats <data_files>");
    exit(0);
    return;
  }
  
  static void displayHelpMessageAndExit(Exception e) {
    writeln(e.msg);
    displayHelpMessageAndExit();
  }
  
  void run() {
    auto M = nrHaplotypes;
    auto nrSubpopPairs = M * (M - 1) / 2;
  
    auto wattersonFactor = TimeIntervals.computeWattersonFactor(M);
    
    size_t calledSites;
    size_t segSites;
    auto pairwiseHets = new double[nrSubpopPairs];
    pairwiseHets[] = 0.0;
    size_t inconsistentInheritance;
    size_t ambiguousPhase;
    auto dummyBoundaries = TimeIntervals.standardTotalBranchlengthIntervals(1, M);
    foreach(filename; fileNames) {
      stderr.writeln("processing file ", filename);
      auto f = File(filename, "r");
      foreach(line; f.byLine) {
        auto fields = line.strip().split();
        calledSites += to!size_t(fields[2]);
        if(canFind(fields[3], '?')) {
          inconsistentInheritance += 1;
          calledSites -= 1;
        }
        else {
          auto alleleStrings = fields[3].split(",");
          bool isSegSite;
          auto nrObs = alleleStrings.length;
          if(nrObs > 1)
            ambiguousPhase += 1;
          foreach(alleleString; alleleStrings) {
            auto pairIndex = 0;
            foreach(i; 0 .. M - 1) {
              foreach(j; i + 1 .. M) {
                if(alleleString[i] != alleleString[j]) {
                  isSegSite = true;
                  pairwiseHets[pairIndex] += 1.0 / nrObs;
                }
                pairIndex += 1;
              }
            }
          }
          if(isSegSite)
            segSites += 1;
        }
      }
    }
    
    writeln("called sites\t",calledSites);
    writeln("segregating sites\t", segSites);
    writeln("inconsistent inheritance\t", inconsistentInheritance);
    writeln("sites with ambiguous phase\t", ambiguousPhase);
    auto pairIndex = 0;
    auto avgTheta = 0.0;
    foreach(i; 0 .. M - 1) {
      foreach(j; i + 1 .. M) {
        writefln("pairwise hets(%s,%s)\t%s", i, j, pairwiseHets[pairIndex]);
        avgTheta += pairwiseHets[pairIndex];
        pairIndex += 1;
      }
    }
    writeln("Tajima's theta\t", avgTheta / (nrSubpopPairs * calledSites));
    writeln("Watterson's theta\t", segSites / (calledSites * wattersonFactor));
  }
}