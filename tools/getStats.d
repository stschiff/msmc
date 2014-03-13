/* Copyright (c) 2012,2013 Genome Research Ltd.
 *
 * Author: Stephan Schiffels <stephan.schiffels@sanger.ac.uk>
 *
 * This file is part of msmc.
 * msmc is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
import std.stdio;
import std.getopt;
import std.string;
import std.conv;
import std.array;
import std.algorithm;
import std.c.stdlib;
import std.exception;
import std.range;

string[] fileNames;
size_t nrHaplotypes;
bool skipAmbiguous;

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
  getopt(args, "skipAmbiguous", &skipAmbiguous);
  fileNames = args[1 .. $];
  enforce(fileNames.length > 0, "need at least one input file");
  nrHaplotypes = getNrHaplotypesFromFile(fileNames[0]);
  enforce(nrHaplotypes >= 2, "need at least two haplotypes");
  stderr.writeln("found ", nrHaplotypes, " haplotypes in file");
}

size_t getNrHaplotypesFromFile(string filename) {
  auto file = File(filename, "r");
  auto line = file.readln();
  line = line.strip();
  auto fields = line.strip().split();
  if(fields.length < 4)
    return 2;
  else {
    auto splitted = fields[3].split(",");
    return cast(size_t)splitted[0].length;
  
  }
}

void printHelp() {
  stderr.writeln("getStats <data_files>");
  exit(0);
  return;
}
  
void run() {
  auto M = nrHaplotypes;
  auto nrSubpopPairs = M * (M - 1) / 2;

  auto wattersonFactor = iota(1, M).map!"1.0/a"().reduce!"a+b"();

  size_t calledSites;
  size_t segSites;
  auto pairwiseHets = new double[nrSubpopPairs];
  pairwiseHets[] = 0.0;
  size_t ambiguousPhase;
  foreach(filename; fileNames) {
    stderr.writeln("processing file ", filename);
    auto f = File(filename, "r");
    foreach(line; f.byLine) {
      auto fields = line.strip().split();
      if(canFind(fields[3], '?'))
        continue;
      calledSites += to!size_t(fields[2]);
      auto alleleStrings = fields[3].split(",");
      auto nrObs = alleleStrings.length;
      if(nrObs > 1)
        ambiguousPhase += 1;
      if(nrObs > 1 && skipAmbiguous)
        continue;
      bool isSegSite;
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
      if(isSegSite) {
        segSites += 1;
      }
    }
  }
  
  writeln("called sites\t",calledSites);
  writeln("segregating sites\t", segSites);
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