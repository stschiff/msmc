import std.exception;
import std.json;
import std.conv;
import std.file;
import std.stdio;
import std.string;
import triple_index_marginal;
import utils;
import time_intervals;
import emission_rate;
import transition_rate;
import coalescence_rate;

class MSMCmodel {
  const EmissionRate emissionRate;
  const TransitionRate transitionRate;
  const MarginalTripleIndex marginalIndex;
  const TimeIntervals timeIntervals;
  const TimeIntervals tTotIntervals;
  const PiecewiseConstantCoalescenceRate coal;

  static MSMCmodel fromFile(string filename) {
    return MSMCmodel.fromJSON(parseJSON(readText(filename)));
  }
  
  static MSMCmodel overrideDemographies(MSMCmodel msmc, string[] filenames) {
    enforce(filenames.length == msmc.nrSubpopulations, "to fix the demographies of all populations, I need as many msmc results as there are populations");
    double newMutationRate, newRecombinationRate;
    auto newLambdaVec = new double[msmc.nrMarginals];
    newLambdaVec[] = 1.0;
    foreach(i; 0 .. msmc.nrSubpopulations) {
      auto json = parseJSON(readText(filenames[i]));
      auto demographySettings = MSMCmodel.fromJSON(json["results"].array[$ - 1]["updatedParams"]);
      enforce(demographySettings.nrSubpopulations == 1, "given demography files must be for exactly one population");
      if(i == 0) {
        newMutationRate = demographySettings.mutationRate;
        newRecombinationRate = demographySettings.recombinationRate;
      }
      auto otherN0overN0 = demographySettings.mutationRate / newMutationRate;
      double[2][] functionValues;
      foreach(j; 0 .. demographySettings.nrTimeIntervals) {
        auto t = demographySettings.timeIntervals.leftBoundary(j) * otherN0overN0;
        auto lambda = demographySettings.lambdaVec[j] / otherN0overN0;
        functionValues ~= [t, lambda];
      }
      auto lambdaVecForSubpop = convertPiecewiseFunctions(msmc.timeIntervals, functionValues);
      
      foreach(au; 0 .. msmc.nrMarginals) {
        auto aij = msmc.marginalIndex.getIndexFromMarginalIndex(au);
        auto triple = msmc.marginalIndex.getTripleFromIndex(aij);
        auto p1 = msmc.subpopLabels[triple.ind1];
        auto p2 = msmc.subpopLabels[triple.ind2];
        if(p1 == i && p2 == i) {
          newLambdaVec[au] = lambdaVecForSubpop[triple.time];
        }
      }
    }
    return new MSMCmodel(newMutationRate, newRecombinationRate, msmc.subpopLabels, newLambdaVec, msmc.timeIntervals.boundaries[0 .. $ - 1], msmc.nrTtotIntervals);
  }
  
  static MSMCmodel fromJSON(JSONValue json) {
    stderr.writeln("creating MSMCmodel from JSON:");
    auto mutationRate = utils.fromJSON!double(json["mutationRate"]);
    auto recombinationRate = utils.fromJSON!double(json["recombinationRate"]);
    auto subpopLabels = utils.fromJSON!(size_t[])(json["subpopLabels"]);
    auto lambdaVec = utils.fromJSON!(double[])(json["lambdaVec"]);
    auto timeBoundaries = utils.fromJSON!(double[])(json["timeIntervals"]);
    auto nrTtotIntervals = utils.fromJSON!size_t(json["nrTtotSegments"]);
    return new MSMCmodel(mutationRate, recombinationRate, subpopLabels, lambdaVec, timeBoundaries, nrTtotIntervals);
  }
  
  this(double mutationRate, double recombinationRate, in size_t[] subpopLabels, in double[] lambdaVec, in double[] timeBoundaries, size_t nrTtotIntervals) {
    auto nrHaplotypes = cast(size_t)subpopLabels.length;
    emissionRate = new EmissionRate(nrHaplotypes, mutationRate);
    timeIntervals = new TimeIntervals(timeBoundaries ~ [double.infinity]);
    tTotIntervals = TimeIntervals.standardTotalBranchlengthIntervals(nrTtotIntervals, nrHaplotypes);
    marginalIndex = new MarginalTripleIndex(nrTimeIntervals, subpopLabels);
    coal = new PiecewiseConstantCoalescenceRate(marginalIndex, lambdaVec);
    transitionRate = new TransitionRate(marginalIndex, coal, timeIntervals, recombinationRate);
  }

  this(double mutationRate, double recombinationRate, size_t[] subpopLabels, double[] lambdaVec,
       size_t nrTimeIntervals, size_t nrTtotIntervals)
  {
    auto standardIntervals = TimeIntervals.standardIntervals(nrTimeIntervals, cast(size_t)subpopLabels.length);
    this(mutationRate, recombinationRate, subpopLabels, lambdaVec, standardIntervals.boundaries[0 .. $ - 1], 
         nrTtotIntervals);
  }
  
  override string toString() const {
    return format("<MSMCmodel: mutationRate=%s, recombinationRate=%s, subpopLabels=%s, lambdaVec=%s, nrTimeIntervals=%s, nrTtotIntervals=%s", mutationRate, recombinationRate, subpopLabels, lambdaVec, nrTimeIntervals, nrTtotIntervals);
  }
  
  static MSMCmodel withTrivialLambda(double mutationRate, double recombinationRate, size_t[] subpopLabels, size_t nrTimeIntervals, size_t nrTtotIntervals) {
    auto marginalIndex = new MarginalTripleIndex(nrTimeIntervals, subpopLabels);
    auto lambdaVec = new double[marginalIndex.nrMarginals];
    lambdaVec[] = 1.0;
    return new MSMCmodel(mutationRate, recombinationRate, subpopLabels, lambdaVec, nrTimeIntervals, nrTtotIntervals);
  }
  
  double emissionProb(string alleles, size_t aij, size_t tTotIndex) const {
    auto triple = marginalIndex.getTripleFromIndex(aij);
    auto type = emissionRate.emissionType(alleles, triple.ind1, triple.ind2);
    auto time = timeIntervals.meanTimeWithLambda(triple.time, coal.getTotalMarginalLambda(triple.time));
    auto tTot = tTotIntervals.meanTime(tTotIndex);
    return emissionRate.emissionProb(type, time, tTot);
  }
  
  double emissionProbHom(size_t time_index, size_t ttotIndex) const {
    auto time = timeIntervals.meanTimeWithLambda(time_index, coal.getTotalMarginalLambda(time_index));
    auto tTot = tTotIntervals.meanTime(ttotIndex);
    return emissionRate.emissionProb(EmissionRate.Observation_t.NoMut, time, tTot);
  }
  
  @property size_t nrHaplotypes() const {
    return cast(size_t)subpopLabels.length;
  }
  
  @property size_t nrMarginals() const {
    return marginalIndex.nrMarginals;
  }
  
  @property size_t nrStates() const {
    return marginalIndex.nrStates;
  }
  
  @property double mutationRate() const {
    return emissionRate.mu;
  }
  
  @property double recombinationRate() const {
    return transitionRate.rho;
  }
  
  @property double[] lambdaVec() const {
    return coal.lambdaVec.dup;
  }
  
  @property size_t[] subpopLabels() const {
    return marginalIndex.subpopLabels.dup;
  }
  
  @property size_t nrSubpopulations() const {
    return marginalIndex.nrSubpopulations;
  }
  
  @property size_t nrTimeIntervals() const {
    return timeIntervals.nrIntervals;
  }
  
  @property size_t nrTtotIntervals() const {
    return tTotIntervals.nrIntervals;
  }
  
  
  JSONValue toJSON() const {
    auto json = JSONValue();
    json.type = JSON_TYPE.OBJECT;
    json.object["mutationRate"] = utils.makeJSON(mutationRate);
    json.object["recombinationRate"] = utils.makeJSON(recombinationRate);
    json.object["timeIntervals"] = utils.makeJSON(timeIntervals.boundaries[0..$ - 1]);
    json.object["nrTtotSegments"] = utils.makeJSON(nrTtotIntervals);
    json.object["lambdaVec"] = utils.makeJSON(lambdaVec);
    json.object["subpopLabels"] = utils.makeJSON(subpopLabels);
    return json;
  }
}

version(unittest) {
  MSMCmodel makeTestParams() {
    auto lambdaVec = [1.2, 3.4, 5.6];
    auto subpopLabels = [0UL,0,0];
    return new MSMCmodel(0.01, 0.001, subpopLabels, lambdaVec, 3, 8);
  }
}

version(unittest) {
  import std.stdio;
  import std.math;
}

unittest {
  writeln("testing ModelParams.toJSON");
  auto p = makeTestParams();
  auto json = p.toJSON();
  assert(approxEqual(json["mutationRate"].floating, 0.01, 1.0e-8, 0.0));
  assert(approxEqual(json["recombinationRate"].floating, 0.001, 1.0e-8, 0.0));
  assert(json["timeIntervals"][0].floating == 0.0);
  assert(json["timeIntervals"][1].floating > 0.0);
  assert(json["timeIntervals"][2].floating > json["timeIntervals"][1].floating);
  assert(json["nrTtotSegments"].integer == 8);
  auto lambdaVec = [1.2, 3.4, 5.6];
  auto subpopLabels = [0U,0,0];
  foreach(i; 0 .. 3)
    assert(approxEqual(json["lambdaVec"][i].floating, lambdaVec[i], 1.0e-12, 0.0));
  foreach(i; 0 .. 3)
    assert(json["subpopLabels"][i].integer == subpopLabels[i]);
}

unittest {
  writeln("testing ModelParams.fromJSON");
  auto p = makeTestParams();
  auto json = p.toJSON();
  auto pNew = MSMCmodel.fromJSON(json);
  assert(pNew.lambdaVec == p.lambdaVec);
  assert(pNew.mutationRate == p.mutationRate);
  assert(pNew.recombinationRate == p.recombinationRate);
  assert(pNew.timeIntervals.boundaries == p.timeIntervals.boundaries);
  assert(pNew.nrTtotIntervals == p.nrTtotIntervals);
  assert(pNew.subpopLabels == p.subpopLabels);
}