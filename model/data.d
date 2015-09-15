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

module model.data; 
import std.stdio;
import std.string;
import std.conv;
import std.algorithm;
import std.math;
import std.c.stdlib;
import std.regex : match, regex, ctRegex;
import std.exception;
import std.range;
import model.time_intervals;

class SegSite_t {
  size_t pos; // rightMost position in the given segment
  size_t[] obs;
  // use vector to allow for alternatives (for example after ambiguous phasing)
  // for missing data: [0]
  
  size_t i_Ttot;
  
  this(size_t pos, in size_t[] obs, size_t i_Ttot) {
    this.pos = pos;
    this.obs = obs.dup;
    this.i_Ttot = i_Ttot;
  }
  
  this(size_t pos, size_t obs, size_t i_Ttot) {
    this.pos = pos;
    this.obs = [obs];
    this.i_Ttot = i_Ttot;
  }
  
  @property SegSite_t dup() const {
    return new SegSite_t(pos, obs.dup, i_Ttot);
  }
  
  override string toString() const {
    return text("Segsite(", pos, ", ", obs, ", ", i_Ttot, ")");
  }
}

string[] canonicalAlleleOrder(size_t M) {
  // we group observations into pairs of strings of two alleles {0,1} which are synonymous with respect to exchanging 0 and 1. For example, we group together the pair [000, 111] or [101, 010]. We denote each group by the version of the string which begins with 0.
  
  string[] allele_order;
  assert(M >= 2);
  auto formatStr = format("%%0%db", M);
  foreach(i; 0 .. 2 ^^ M) {
    allele_order ~= format(formatStr, i);
  }
  return allele_order;
}

unittest {
  writeln("test canonicalAlleleOrder");
  assert(canonicalAlleleOrder(2) == ["00", "01", "10", "11"]);
  assert(canonicalAlleleOrder(3) == ["000", "001", "010", "011", "100", "101", "110", "111"]);
}

string invertAllele(string allele) {
  auto newA = new char[allele.length];
  foreach(i, a; allele) {
    assert(a == '0' || a == '1');
    newA[i] = a == '0' ? '1' : '0';
  }
  return newA.idup;
}

unittest {
  writeln("test invertAllele");
  assert(invertAllele("0011") == "1100");
  assert(invertAllele("1111") == "0000");
  assert(invertAllele("0010") == "1101");
  assert(invertAllele("01") == "10");
}


// this function converts any allele string to a sequence of 0's and 1's
string normalizeAlleleString(string alleles) {
  char firstAllele = alleles[0];
  char[] ret;
  foreach(allele; alleles) {
    ret ~= allele == firstAllele ? '0' : '1';
  }
  return ret.idup;
}

unittest {
  writeln("test normalizeAlleleString");
  assert(normalizeAlleleString("AACC") == "0011");
  assert(normalizeAlleleString("1101") == "0010");
  assert(normalizeAlleleString("ACTG") == "0111");
  assert(normalizeAlleleString("GGGG") == "0000");
  assert(normalizeAlleleString("GG") == "00");
  assert(normalizeAlleleString("TC") == "01");
}

// void checkDataLine(const char[] line) {
//   auto r = regex(r"^\w+\s\d+\s\d+(\s[ACTG01\?,]+){0,1}$");
//   enforce(match(line, r));
// }
//
// unittest {
//   assertThrown(checkDataLine("1 20 5 AACC,AACA 2.44"));
//   assertNotThrown(checkDataLine("1 20 5 AACC"));
//   assertNotThrown(checkDataLine("4 5 2"));
//   assertNotThrown(checkDataLine("1 10 5 ACC"));
//   assertThrown(checkDataLine("1 20 5 AGGSSXX"));
// }

size_t getNrHaplotypesFromFile(string filename) {
  auto f = File(filename, "r");
  scope(exit) f.close();
  auto line = f.readln();
  line = line.strip();
  // checkDataLine(line);
  auto fields = line.strip().split();
  enforce(fields.length > 2, "illegal input file");
  if(fields.length < 4)
    return 2;
  else {
    auto splitted = fields[3].split(",");
    return cast(size_t)splitted[0].length;
  }
}

unittest {
  auto tmp = File("/tmp/nrHaplotypesTest.txt", "w");
  tmp.writeln("1 10 5 ACC,CCA");
  tmp.close();
  assert(getNrHaplotypesFromFile("/tmp/nrHaplotypesTest.txt") == 3);
  tmp = File("/tmp/nrHaplotypesTest.txt", "w");
  tmp.writeln("1 10 5");
  tmp.close();
  assert(getNrHaplotypesFromFile("/tmp/nrHaplotypesTest.txt") == 2);
}

SegSite_t[] readSegSites(string filename, bool directedEmissions, size_t[] indices, bool skipAmbiguous) {
  // format: chr, position, nr_calledSites, [alleles]
  // if no alleles are given, assume M=2 and "01"
  // alleles can be given as comma-separated list of alternative alleles
  
  SegSite_t[] ret;

  size_t M;
  if(indices.length == 0) {
    M = getNrHaplotypesFromFile(filename);
    indices = iota(M).array();
  }
  else {
    M = indices.length;
  }
  int obsMap[string];
  auto allele_order = canonicalAlleleOrder(M);
  auto index = 1; // index=0 indicates missing data 
  foreach(allele; allele_order) {
    obsMap[allele] = index;
    ++index;
  }
  
  auto f = File(filename, "r");
  long lastPos = -1;
  auto chrom = "";
  foreach(line; f.byLine()) {
    // checkDataLine(line.strip());
    auto fields = line.strip().split();
    if(chrom == "")
        chrom = fields[0].idup;
    else
        enforce(chrom == fields[0], "chromosomes must all be the same within one file (sorry)");
    auto pos = to!size_t(fields[1]);
    auto nrCalledSites = to!size_t(fields[2]);
    enforce(nrCalledSites > 0, "nr of called sites (3rd column in input file) must be always > 0");
    if(lastPos == -1) {
      lastPos = pos - nrCalledSites;
    }
    
    enforce(nrCalledSites <= pos - lastPos);
    enforce(nrCalledSites > 0);
    
    if(fields.length > 2) {
      // checking whether we have any "N" or "?" in the data, which would mark it as missing data.
      auto is_missing = false;
      auto raw_allele_strings = split(fields[3], ",");
      if(!is_missing) {
        foreach(raw_allele_string; raw_allele_strings) {
          foreach(i; indices) {
            if(i >= raw_allele_string.length) {
              stderr.writeln("Haplotype index exceeds number of haplotypes in datafile");
              exit(0);
            }
            if(!canFind("ACTG01", raw_allele_string[i])) {
              is_missing = true;
              break;
            }
          }
          // if(indices.any!(a => !canFind("ACTG01", raw_allele_string[a]))) {
          //   is_missing = true;
          //   break;
          // }
        }
      }
      if(is_missing) {
        if(nrCalledSites < pos - lastPos) { // missing data
          ret ~= new SegSite_t(pos - nrCalledSites, 0, 0);
        }
        if(nrCalledSites > 1)
          ret ~= new SegSite_t(pos - 1, 1, 0);
        ret ~= new SegSite_t(pos, 0, 0);
        lastPos = pos;
      }
      else {
        size_t[] allele_indices;
        foreach(allele_string; split(fields[3], ",")) {
          char[] selected_allele_string;
          foreach(i; indices)
            selected_allele_string ~= allele_string[i];
          enforce(selected_allele_string.length == M);
          auto normalized = directedEmissions ? selected_allele_string : 
                            normalizeAlleleString(selected_allele_string.idup);
          allele_indices ~= obsMap[normalized];
        }
        if(nrCalledSites < pos - lastPos) { // missing data
          ret ~= new SegSite_t(pos - nrCalledSites, 0, 0);
        }
        allele_indices = allele_indices.uniq().array();
        if(skipAmbiguous && allele_indices.length > 1)
          ret ~= new SegSite_t(pos, 0, 0);
        else
          ret ~= new SegSite_t(pos, allele_indices, 0);
        lastPos = pos;
      }
    }
    else {
      if(nrCalledSites < pos - lastPos) { // missing data
        ret ~= new SegSite_t(pos - nrCalledSites, 0, 0);
      }
      ret ~= new SegSite_t(pos, 2, 0); // [2] means heterozygous
      lastPos = pos;
    }
  }
  
  foreach(i; 1 .. ret.length) {
    assert(ret[i].pos > ret[i - 1].pos, text([i, ret[i].pos, ret[i - 1].pos]));
  }
  
  enforce(ret.length > 0, "found empty input file: " ~ filename);
  return ret;
}

unittest {
  writeln("test readSegSites");
  auto tmp_file = File("/tmp/msmc_data_unittest.tmp", "w");
  tmp_file.writeln("1 1000000 42 AACC");
  tmp_file.writeln("1 1000004 2 ACCG");
  tmp_file.writeln("1 1000008 3 ACC?,ATTA");
  tmp_file.writeln("1 1000012 4 ACCG,TTGA");
  tmp_file.close();

  auto segsites = readSegSites("/tmp/msmc_data_unittest.tmp", false, [], false);
  assert(segsites[0].pos == 1000000 && segsites[0].obs == [4]);
  assert(segsites[1].pos == 1000002 && segsites[1].obs == [0]);
  assert(segsites[3].pos == 1000005 && segsites[3].obs == [0]);
  assert(segsites[4].pos == 1000007 && segsites[4].obs == [1]);
  assert(segsites[5].pos == 1000008 && segsites[5].obs == [0]);
  assert(segsites[6].pos == 1000012 && segsites[6].obs == [8, 4]);
  
}

unittest {
  writeln("test readSegSites for M=2");
  auto tmp_file = File("/tmp/msmc_data_unittest.tmp", "w");
  tmp_file.writeln("1 1000000 42 AC");
  tmp_file.writeln("1 1000004 2 CC");
  tmp_file.writeln("1 1000008 3 C?,AT");
  tmp_file.writeln("1 1000012 4 AA,TT");
  tmp_file.close();

  auto segsites = readSegSites("/tmp/msmc_data_unittest.tmp", false, [], false);
  assert(segsites[0].pos == 1000000 && segsites[0].obs == [2]);
  assert(segsites[1].pos == 1000002 && segsites[1].obs == [0]);
  assert(segsites[3].pos == 1000005 && segsites[3].obs == [0]);
  assert(segsites[4].pos == 1000007 && segsites[4].obs == [1]);
  assert(segsites[5].pos == 1000008 && segsites[5].obs == [0]);
  assert(segsites[6].pos == 1000012 && segsites[6].obs == [1]);
  
}

double getTheta(in SegSite_t[][] data, size_t nrHaplotypes) {
  size_t nr_segsites; 
  size_t called_sites;
  foreach(d; data) {
    size_t lastPos = 0;
    foreach(dd; d) {
      if(dd.obs[0] > 0) {
        if(lastPos > 0)
          called_sites += dd.pos - lastPos;
        if(dd.obs.any!(o => o > 1))
          nr_segsites += 1;
      }
      lastPos = dd.pos;
    }
  }
  auto watterson = iota(1, nrHaplotypes).map!"1.0 / a"().reduce!"a+b"();
  return cast(double)nr_segsites / cast(double)called_sites / watterson;
}