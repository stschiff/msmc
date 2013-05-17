import std.stdio;
import std.exception;
import triple_index;
import triple_index_marginal;

class PiecewiseConstantCoalescenceRate {
  const double[] lambdaVec;
private:
  size_t nrTimeSegments, nrHaplotypes;
  double[][][] lambda;
  double[] totalMarginalLambda; 
  double[] avgLambda;
  double[][] semiMarginalLambda;
    
public:
  this(in MarginalTripleIndex marginalIndex, in double[] lambdaVec) {
    enforce(lambdaVec.length == marginalIndex.nrMarginals);
    this.lambdaVec = lambdaVec;
    foreach(l; lambdaVec)
      enforce(l > 0.0, "need positive lambdaVec");
    nrTimeSegments = marginalIndex.nrTimeIntervals;
    nrHaplotypes = marginalIndex.nrIndividuals;
    lambda = new double[][][](nrTimeSegments, nrHaplotypes, nrHaplotypes);
    totalMarginalLambda = new double[nrTimeSegments];
    avgLambda = new double[nrTimeSegments];
    semiMarginalLambda = new double[][](nrTimeSegments, nrHaplotypes);
    
    foreach(t; 0 .. nrTimeSegments) {
      foreach(i; 0 .. nrHaplotypes) foreach(j; 0 .. nrHaplotypes) {
        auto triple = Triple(t, i, j);
        auto mIndex = marginalIndex.getMarginalIndexFromGeneralTriple(triple);
        lambda[t][i][j] = lambdaVec[mIndex];
      }
      totalMarginalLambda[t] = computeTotalMarginalLambda(t);
      avgLambda[t] = computeAvgLambda(t);
      foreach(i; 0 .. nrHaplotypes)
        semiMarginalLambda[t][i] = computeSemiMarginalLambda(t, i);
    }
  }
  
  double getLambda(size_t time, size_t ind1, size_t ind2) const {
    return lambda[time][ind1][ind2];
  }
  
  double getTotalMarginalLambda(size_t time) const {
    return totalMarginalLambda[time];
  }
  
  double getAvgLambda(size_t time) const {
    return avgLambda[time];
  }
  
  double getSemiMarginalLambda(size_t time, size_t branch) const {
    return semiMarginalLambda[time][branch];
  }
  
  private double computeTotalMarginalLambda(size_t time) const {
    auto sum = 0.0;
    foreach(ind1; 0 .. nrHaplotypes - 1) foreach(ind2; ind1 + 1 .. nrHaplotypes)
      sum += getLambda(time, ind1, ind2);
    return sum;
  }
  
  private double computeAvgLambda(size_t time) const {
    auto mOver2 = nrHaplotypes * (nrHaplotypes - 1) / 2;
    return getTotalMarginalLambda(time) / mOver2;
  }
  
  private double computeSemiMarginalLambda(size_t time, size_t branch) const {
    auto sum = 0.0;
    foreach(branch2; 0 .. nrHaplotypes) {
      sum += getLambda(time, branch, branch2);
    }
    return sum;
  }
}

version(unittest) {
  PiecewiseConstantCoalescenceRate makeTestCoal() {
    auto subPopLabels = [0UL, 1, 0, 1];
    auto marginalIndex = new MarginalTripleIndex(2UL, subPopLabels);
    auto lambdaVec = [1.0, 1.5, 2.0, 3.0, 5.0, 3.0];
    return new PiecewiseConstantCoalescenceRate(marginalIndex, lambdaVec);
  }
  PiecewiseConstantCoalescenceRate makeTestCoalM2() {
    auto subPopLabels = [0UL, 0];
    auto marginalIndex = new MarginalTripleIndex(2U, subPopLabels);
    auto lambdaVec = [1.0, 1.5];
    return new PiecewiseConstantCoalescenceRate(marginalIndex, lambdaVec);
  }
}

unittest {
  writeln("test coal.getLambda");
  auto coal = makeTestCoal();
  assert(coal.getLambda(0, 0, 1) == 1.5);
  assert(coal.getLambda(1, 1, 2) == 5.0);
  assert(coal.getLambda(1, 0, 2) == 3.0);
  coal = makeTestCoalM2();
  assert(coal.getLambda(0, 0, 0) == 1.0);
  assert(coal.getLambda(0, 0, 1) == 1.0);
  assert(coal.getLambda(0, 1, 0) == 1.0);
  assert(coal.getLambda(0, 1, 1) == 1.0);
  assert(coal.getLambda(1, 0, 0) == 1.5);
  assert(coal.getLambda(1, 0, 1) == 1.5);
  assert(coal.getLambda(1, 1, 0) == 1.5);
  assert(coal.getLambda(1, 1, 1) == 1.5);
}

unittest {
  writeln("test coal.getTotalMarginalLambda");
  auto coal = makeTestCoal();
  assert(coal.getTotalMarginalLambda(0) == 9.0);
  coal = makeTestCoalM2();
  assert(coal.getTotalMarginalLambda(0) == 1.0);
  assert(coal.getTotalMarginalLambda(1) == 1.5);
}

unittest {
  writeln("test coal.getAvgLambda");
  auto coal = makeTestCoal();
  assert(coal.getAvgLambda(0) == 1.5);
  coal = makeTestCoalM2();
  assert(coal.getAvgLambda(0) == 1.0);
  assert(coal.getAvgLambda(1) == 1.5);
}

unittest {
  writeln("test coal.getSemiMarginalLambda");
  auto coal = makeTestCoal();
  assert(coal.getSemiMarginalLambda(0, 0) == 5.0);
  assert(coal.getSemiMarginalLambda(0, 1) == 7.0);
  coal = makeTestCoalM2();
  assert(coal.getSemiMarginalLambda(0, 0) == 2.0);
  assert(coal.getSemiMarginalLambda(0, 1) == 2.0);
  assert(coal.getSemiMarginalLambda(1, 0) == 3.0);
  assert(coal.getSemiMarginalLambda(1, 1) == 3.0);
  
}
