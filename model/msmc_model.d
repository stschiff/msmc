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
 
module model.msmc_model;
import std.exception;
import std.json;
import std.conv;
import std.file;
import std.stdio;
import std.string;
import model.triple_index_marginal;
import model.time_intervals;
import model.emission_rate;
import model.transition_rate;
import model.coalescence_rate;

class MSMCmodel {
  const EmissionRate emissionRate;
  const TransitionRate transitionRate;
  const MarginalTripleIndex marginalIndex;
  const TimeIntervals timeIntervals;
  const TimeIntervals tTotIntervals;
  const PiecewiseConstantCoalescenceRate coal;

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
    // auto time = timeIntervals.meanTimeWithLambda(triple.time, coal.getTotalMarginalLambda(triple.time));
    auto time = timeIntervals.meanTimeWithLambda(triple.time, nrHaplotypes);
    auto tTot = tTotIntervals.meanTime(tTotIndex, nrHaplotypes);
    return emissionRate.emissionProb(type, time, tTot);
  }
  
  double emissionProbHom(size_t time_index, size_t ttotIndex) const {
    // auto time = timeIntervals.meanTimeWithLambda(time_index, coal.getTotalMarginalLambda(time_index));
    auto time = timeIntervals.meanTimeWithLambda(time_index, nrHaplotypes);
    auto tTot = tTotIntervals.meanTime(ttotIndex, nrHaplotypes);
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

// unittest {
//   writeln("testing ModelParams.toJSON");
//   auto p = makeTestParams();
//   auto json = p.toJSON();
//   assert(approxEqual(json["mutationRate"].floating, 0.01, 1.0e-8, 0.0));
//   assert(approxEqual(json["recombinationRate"].floating, 0.001, 1.0e-8, 0.0));
//   assert(json["timeIntervals"][0].floating == 0.0);
//   assert(json["timeIntervals"][1].floating > 0.0);
//   assert(json["timeIntervals"][2].floating > json["timeIntervals"][1].floating);
//   assert(json["nrTtotSegments"].integer == 8);
//   auto lambdaVec = [1.2, 3.4, 5.6];
//   auto subpopLabels = [0U,0,0];
//   foreach(i; 0 .. 3)
//     assert(approxEqual(json["lambdaVec"][i].floating, lambdaVec[i], 1.0e-12, 0.0));
//   foreach(i; 0 .. 3)
//     assert(json["subpopLabels"][i].integer == subpopLabels[i]);
// }
// 
// unittest {
//   writeln("testing ModelParams.fromJSON");
//   auto p = makeTestParams();
//   auto json = p.toJSON();
//   auto pNew = MSMCmodel.fromJSON(json);
//   assert(pNew.lambdaVec == p.lambdaVec);
//   assert(pNew.mutationRate == p.mutationRate);
//   assert(pNew.recombinationRate == p.recombinationRate);
//   assert(pNew.timeIntervals.boundaries == p.timeIntervals.boundaries);
//   assert(pNew.nrTtotIntervals == p.nrTtotIntervals);
//   assert(pNew.subpopLabels == p.subpopLabels);
// }