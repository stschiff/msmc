import std.algorithm;
import std.range;
import std.conv;
import std.stdio;
import std.exception;
import triple_index;

unittest {
  writeln("test MarginalIndex.getMarginalIndex");
  auto mI = new MarginalTripleIndex(4, [0, 1, 0, 1]);
  assert(mI.getMarginalIndexFromGeneralTriple(Triple(0, 0, 1)) == 1);
  assert(mI.getMarginalIndexFromGeneralTriple(Triple(0, 0, 2)) == 0);
  assert(mI.getMarginalIndexFromGeneralTriple(Triple(0, 2, 3)) == 1);
  assert(mI.getMarginalIndexFromGeneralTriple(Triple(1, 0, 2)) == 3);
  assert(mI.getMarginalIndexFromGeneralTriple(Triple(3, 2, 3)) == 10);
}

unittest {
  writeln("test MarginalIndex.getTriple");
  auto mI = new MarginalTripleIndex(4, [0, 1, 0, 1]);
  assert(mI.getDegeneracyForMarginalIndex(1) == 4);
  auto index = mI.getIndexFromMarginalIndex(1, 2);
  auto triple = mI.getTripleFromIndex(index);
  assert(triple == Triple(0, 1, 2));
}

class MarginalTripleIndex : TripleIndex {
  private size_t nrSubpopulations_;
  size_t[] subpopLabels;
  
  private size_t[][] marginalIndexToIndicesMap;
  private size_t[] indexToMarginalIndexMap;
  private size_t[][][] subpopulationTripleToMarginalIndexMap;
  
  this(size_t nrTimeIntervals, in size_t[] subpopLabels)
  {
    enforce(subpopLabels.length >= 2, "need at least two haplotypes");
    enforce(nrTimeIntervals > 0);
    auto nrIndividuals = cast(size_t)subpopLabels.length;
    super(nrTimeIntervals, nrIndividuals);
    this.subpopLabels = subpopLabels.dup;
    nrSubpopulations_ = computeNrSubpops(subpopLabels);
    
    computeMaps();
  }

  static size_t computeNrSubpops(in size_t[] subpopLabels) {
    auto labels = subpopLabels.dup;
    auto uniqueLabels = uniq(sort(labels));
    auto nr = cast(size_t)walkLength(uniqueLabels);
    auto maxLabel = minCount!"a>b"(labels)[0];
    enforce(nr == maxLabel + 1, text(nr, " ", maxLabel));
    return nr;
  }
  
  private void computeMaps() {
    computeSubpopulationTripleMap();
    computeIndexToIndexMaps();
  }
  
  private void computeSubpopulationTripleMap() {
    subpopulationTripleToMarginalIndexMap = new size_t[][][](nrTimeIntervals, nrSubpopulations_, nrSubpopulations_);
    auto marginalIndex = 0U;
    foreach(time; 0 .. nrTimeIntervals) {
      foreach(subpop1; 0 .. nrSubpopulations_) {
        foreach(subpop2; subpop1 .. nrSubpopulations_) { 
          subpopulationTripleToMarginalIndexMap[time][subpop1][subpop2] = marginalIndex;
          subpopulationTripleToMarginalIndexMap[time][subpop2][subpop1] = marginalIndex;
          ++marginalIndex;
        }
      }
    }
  }
  
  private void computeIndexToIndexMaps() {
    marginalIndexToIndicesMap = new size_t[][nrMarginals()];
    indexToMarginalIndexMap = new size_t[nrStates()];
    foreach(index; 0 .. nrStates()) {
      auto triple = getTripleFromIndex(index);
      auto subpop1 = subpopLabels[triple.ind1];
      auto subpop2 = subpopLabels[triple.ind2];
      auto marginalIndex = subpopulationTripleToMarginalIndexMap[triple.time][subpop1][subpop2];
      indexToMarginalIndexMap[index] = marginalIndex;
      marginalIndexToIndicesMap[marginalIndex] ~= index;
    }
  }
  
  @property size_t nrMarginals() const {
    return nrTimeIntervals * nrSubpopulations_ * (nrSubpopulations_ + 1) / 2;
  }
  
  @property size_t nrSubpopulations() const {
    return nrSubpopulations_;
  }
  
  size_t getMarginalIndexFromIndex(size_t index) const {
    return indexToMarginalIndexMap[index];
  }
  
  size_t getMarginalIndexFromGeneralTriple(Triple triple) const
    in {
      assert(triple.time < nrTimeIntervals &&
          triple.ind1 < nrIndividuals && triple.ind2 < nrIndividuals, text(triple));
    }
  body {
    auto subpop1 = cast()subpopLabels[triple.ind1];
    auto subpop2 = cast()subpopLabels[triple.ind2];
    return subpopulationTripleToMarginalIndexMap[triple.time][subpop1][subpop2];
  }
  
  size_t getIndexFromMarginalIndex(size_t marginalIndex, size_t degeneracyIndex=0) const {
    return marginalIndexToIndicesMap[marginalIndex][degeneracyIndex];
  }
  
  size_t getDegeneracyForMarginalIndex(size_t marginalIndex) const {
    return cast(size_t)marginalIndexToIndicesMap[marginalIndex].length;
  }
  
}

