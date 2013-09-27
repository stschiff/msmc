import std.stdio;
import std.getopt;
import std.exception;
import std.file;
import std.conv;
import std.string;
import std.algorithm;
import std.array;
import std.range;
import std.typecons;

string[] haplotypeFileNames;
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
  getopt(args, "hapsFile|h", &haplotypeFileNames, "multihetsepIndex|i", &multihetsepIndices, "hapsIndex|j", &hapsIndices);
  nrIndividuals = multihetsepIndices.length;
  enforce(nrIndividuals > 0, "need to specify at least one individual");
  enforce(hapsIndices.length == nrIndividuals && haplotypeFileNames.length == nrIndividuals, "need same number of indices in multihetsep and hapsfile and correct number of hapsFiles");
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
  auto multihetsepData = readMultihetsepFile(multihetsep_filename);
  
  foreach(ind_index; 0 .. nrIndividuals) {
    auto hapsFile = File(haplotypeFileNames[ind_index], "r");
    auto hapsIndex = hapsIndices[ind_index];
    auto multihetsepIndex = multihetsepIndices[ind_index];
    auto lastHapsFilePos = 0UL;
    auto hapsFileRange = hapsFile.byLine;
    char[][] hapsFields;
    foreach(ref data; multihetsepData) {
      auto chr = data[0];
      auto pos = data[1];
      auto nr_called = data[2];
      auto unphased = data[3];
      auto allPhasings = data[4];
      while(lastHapsFilePos < pos && !hapsFileRange.empty) {
        auto hapsFileLine = hapsFileRange.front.dup;
        hapsFileRange.popFront;
        hapsFields = hapsFileLine.strip.split;
        lastHapsFilePos = hapsFields[2].to!size_t();
      }
      if(lastHapsFilePos == pos) {
        auto al = hapsFields[3..5].joiner.array;
        auto gens = hapsFields[5 + 2 * hapsIndex .. 5 + 2 * (hapsIndex + 1)].map!(a => al[a.to!size_t()].to!char())().array;
        auto prunedPhasings = pruneInconsistentPhasings(allPhasings, multihetsepIndex, gens);
        enforce(prunedPhasings.length > 0, text(format("%s\t%s\t%s\t%s\t%s", pos, unphased, allPhasings, al, gens)));
        data[4] = prunedPhasings;
      }
    }
  }
  foreach(data; multihetsepData) {
    writefln("%s\t%s\t%s\t%s", data[0], data[1], data[2], data[4].joiner(",").array());
  }
}

auto readMultihetsepFile(string filename) {
  Tuple!(string, size_t, size_t, char[][], char[][])[] ret;
  auto f = File(multihetsep_filename, "r");
  foreach(line; f.byLine) {
    auto fields = line.strip().split();
    auto chr = fields[0].idup;
    auto pos = fields[1].to!size_t();
    auto calledSites = fields[2].to!size_t();
    auto unphased = fields[3 .. $];
    auto allPhasings = getAllPhasings(unphased);
    ret ~= tuple(chr, pos, calledSites, unphased, allPhasings);
  }
  return ret;
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

char[][] pruneInconsistentPhasings(in char[][] allPhasings, size_t ind, in char[] phasedGenotypes) {
  char[][] ret;
  foreach(alleles; allPhasings) {
    if(alleles[ind * 2] == phasedGenotypes[0] && alleles[ind * 2 + 1] == phasedGenotypes[1])
      ret ~= alleles.dup;
  }
  return ret;
}

unittest {
  auto seqs = ["ACAAAC", "ACAACA", "CAAAAC", "CAAACA"].map!"a.dup"().array();
  auto pruned = pruneInconsistentPhasings(seqs, 0, ['C', 'A']);
  assert(equal(pruned, ["CAAAAC", "CAAACA"]), text(pruned));
  pruned = pruneInconsistentPhasings(seqs, 2, ['A', 'C']);
  assert(equal(pruned, ["ACAAAC", "CAAAAC"]), text(pruned));
}
