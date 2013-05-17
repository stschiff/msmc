import std.stdio;
import std.string;
import std.conv;
import std.algorithm;
import std.math;
import std.c.stdlib;
import std.regex : match, regex, ctRegex;
import std.exception;
import time_intervals;

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
  foreach(i; 0 .. 2 ^^ (M - 1)) {
    allele_order ~= format(formatStr, i);
  }
  return allele_order;
}

unittest {
  writeln("test canonicalAlleleOrder");
  assert(canonicalAlleleOrder(2) == ["00", "01"]);
  assert(canonicalAlleleOrder(3) == ["000", "001", "010", "011"]);
  assert(canonicalAlleleOrder(4) == [
      "0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111"
    ]);
  assert(canonicalAlleleOrder(5) == [
      "00000", "00001", "00010", "00011", "00100", "00101", "00110", "00111",
      "01000", "01001", "01010", "01011", "01100", "01101", "01110", "01111"
    ]);
  assert(canonicalAlleleOrder(6) == [
      "000000", "000001", "000010", "000011", "000100", "000101", "000110", "000111",
      "001000", "001001", "001010", "001011", "001100", "001101", "001110", "001111",
      "010000", "010001", "010010", "010011", "010100", "010101", "010110", "010111",
      "011000", "011001", "011010", "011011", "011100", "011101", "011110", "011111"
    ]);
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

void checkDataLine(const char[] line) {
  auto r = regex(r"^\w+\s\d+\s\d+(\s[ACTG01\?,]+(\s[\d\.]+){0,1}){0,1}$");
  // enum r = ctRegex!(r"^\w+\s\d+\s\d+(\s[ACTG01\?,]+(\s[\d\.]+){0,1}){0,1}$");
  enforce(match(line, r));
}

unittest {
  assertNotThrown(checkDataLine("1 20 5 AACC,AACA 4.56"));
  assertNotThrown(checkDataLine("1 20 5 AACC"));
  assertNotThrown(checkDataLine("4 5 2"));
  assertNotThrown(checkDataLine("1 10 5 ACC 2.45"));
  assertThrown(checkDataLine("1 20 5 AGGSSXX 5.11"));
  assertThrown(checkDataLine("1 20 5 5.11"));
}

size_t getNrHaplotypesFromFile(string filename) {
  scope(exit) file.close();
  auto file = File(filename, "r");
  auto line = file.readln();
  line = line.strip();
  checkDataLine(line);
  auto fields = line.strip().split();
  if(fields.length < 4)
    return 2;
  else {
    auto splitted = fields[3].split(",");
    return cast(size_t)splitted[0].length;
  
  }
}

unittest {
  auto tmp = File("/tmp/nrHaplotypesTest.txt", "w");
  tmp.writeln("1 10 5 ACC,CCA 2.45");
  tmp.close();
  assert(getNrHaplotypesFromFile("/tmp/nrHaplotypesTest.txt") == 3);
  tmp = File("/tmp/nrHaplotypesTest.txt", "w");
  tmp.writeln("1 10 5");
  tmp.close();
  assert(getNrHaplotypesFromFile("/tmp/nrHaplotypesTest.txt") == 2);
}

SegSite_t[] readSegSites(string filename, size_t M, in TimeIntervals ttotIntervals) {
  // format: chr, position, nr_calledSites, [alleles, [ttot]]
  // if no alleles are given, assume M=2 and "01"
  // alleles can be given as comma-separated list of alternative alleles
  
  stderr.writeln("reading data from file: ", filename);
  assert(M >= 2);
  
  SegSite_t[] ret;

  int obsMap[string];
  auto allele_order = canonicalAlleleOrder(M);
  auto index = 1; // index=0 indicates missing data 
  foreach(allele; allele_order) {
    obsMap[allele] = index;
    obsMap[invertAllele(allele)] = index; // symmetrizing states: 1101 is the same as 0010 !
    ++index;
  }
  
  auto f = File(filename, "r");
  long lastPos = -1;
  foreach(line; f.byLine()) {
    // checkDataLine(line.strip());
    auto fields = line.strip().split();
    auto pos = to!size_t(fields[1]);
    auto nrCalledSites = to!size_t(fields[2]);
    if(lastPos == -1) {
      lastPos = pos - nrCalledSites;
    }
    
    // use nan as dummy value for Ttot
    double Ttot; 
    if(fields.length > 4)
      Ttot = to!double(fields[4]);
    else
      Ttot = 0.0;
    auto i_Ttot = ttotIntervals.findIntervalForTime(Ttot);
    
    enforce(nrCalledSites <= pos - lastPos);
    enforce(nrCalledSites > 0);
    
    if(fields.length > 2) {
      // checking whether we have any "N" or "?" in the data, which would mark it as missing data.
      auto is_missing = false;
      foreach(raw_allele_string; split(fields[3], ",")) {
        foreach(pos_; 0 .. M) {
          if(pos_ >= raw_allele_string.length) {
            stderr.writefln("Haplotype index %s exceeds number of haplotypes in datafile", pos_);
            exit(0);
          }
          if(!canFind("ACTG01", raw_allele_string[pos_])) {
            is_missing = true;
            break;
          }
        }
      }
      if(is_missing) {
        if(nrCalledSites < pos - lastPos) { // missing data
          ret ~= new SegSite_t(pos - nrCalledSites, 0, i_Ttot);
        }
        if(nrCalledSites > 1)
          ret ~= new SegSite_t(pos - 1, 1, i_Ttot);
        ret ~= new SegSite_t(pos, 0, i_Ttot);
        lastPos = pos;
      }
      else {
        size_t[] allele_indices;
        foreach(allele_string; split(fields[3], ",")) {
          enforce(allele_string.length == M);
          auto normalized = normalizeAlleleString(allele_string.idup);
          allele_indices ~= obsMap[normalized];
        }
        if(nrCalledSites < pos - lastPos) { // missing data
          ret ~= new SegSite_t(pos - nrCalledSites, 0, i_Ttot);
        }
        ret ~= new SegSite_t(pos, allele_indices, i_Ttot);
        lastPos = pos;
      }
    }
    else {
      if(nrCalledSites < pos - lastPos) { // missing data
        ret ~= new SegSite_t(pos - nrCalledSites, 0, i_Ttot);
      }
      ret ~= new SegSite_t(pos, 2, i_Ttot); // [2] means heterozygous
      lastPos = pos;
    }
    
  }
  
  foreach(i; 1 .. ret.length) {
    assert(ret[i].pos > ret[i - 1].pos, text([i, ret[i].pos, ret[i - 1].pos]));
  }
  
  return ret;
}

unittest {
  writeln("test readSegSites");
  auto tmp_file = File("/tmp/msmc_data_unittest.tmp", "w");
  tmp_file.writeln("1 1000000 42 AACC 2.56");
  tmp_file.writeln("1 1000004 2 ACCG 2.3");
  tmp_file.writeln("1 1000008 3 ACC?,ATTA 4.55");
  tmp_file.writeln("1 1000012 4 ACCG,TTGA 9.123");
  tmp_file.close();

  auto ttotIntervals = TimeIntervals.standardTotalBranchlengthIntervals(1, 4);
  auto segsites = readSegSites("/tmp/msmc_data_unittest.tmp", 4, ttotIntervals);
  foreach(s; segsites)
    assert(s.i_Ttot == 0);
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

  auto ttotIntervals = TimeIntervals.standardTotalBranchlengthIntervals(1, 2);
  auto segsites = readSegSites("/tmp/msmc_data_unittest.tmp", 2, ttotIntervals);
  foreach(s; segsites)
    assert(s.i_Ttot == 0);
  assert(segsites[0].pos == 1000000 && segsites[0].obs == [2]);
  assert(segsites[1].pos == 1000002 && segsites[1].obs == [0]);
  assert(segsites[3].pos == 1000005 && segsites[3].obs == [0]);
  assert(segsites[4].pos == 1000007 && segsites[4].obs == [1]);
  assert(segsites[5].pos == 1000008 && segsites[5].obs == [0]);
  assert(segsites[6].pos == 1000012 && segsites[6].obs == [1, 1]);
  
}