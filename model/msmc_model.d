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

  this(double mutationRate, double recombinationRate, in size_t[] subpopLabels, in double[] lambdaVec, in double[] timeBoundaries, size_t nrTtotIntervals, bool directedEmissions) {
    auto nrHaplotypes = cast(size_t)subpopLabels.length;
    emissionRate = new EmissionRate(nrHaplotypes, mutationRate, directedEmissions);
    timeIntervals = new TimeIntervals(timeBoundaries ~ [double.infinity]);
    tTotIntervals = TimeIntervals.standardTotalBranchlengthIntervals(nrTtotIntervals);
    marginalIndex = new MarginalTripleIndex(nrTimeIntervals, subpopLabels);
    coal = new PiecewiseConstantCoalescenceRate(marginalIndex, lambdaVec);
    transitionRate = new TransitionRate(marginalIndex, coal, timeIntervals, recombinationRate);
  }

  this(double mutationRate, double recombinationRate, size_t[] subpopLabels, double[] lambdaVec,
       size_t nrTimeIntervals, size_t nrTtotIntervals, bool directedEmissions)
  {
    auto standardIntervals = TimeIntervals.standardIntervals(nrTimeIntervals, cast(size_t)subpopLabels.length);
    this(mutationRate, recombinationRate, subpopLabels, lambdaVec, standardIntervals.boundaries[0 .. $ - 1], 
         nrTtotIntervals, directedEmissions);
  }
  
  override string toString() const {
    return format("<MSMCmodel: mutationRate=%s, recombinationRate=%s, subpopLabels=%s, lambdaVec=%s, nrTimeIntervals=%s, nrTtotIntervals=%s", mutationRate, recombinationRate, subpopLabels, lambdaVec, nrTimeIntervals, nrTtotIntervals);
  }
  
  static MSMCmodel withTrivialLambda(double mutationRate, double recombinationRate, size_t[] subpopLabels, size_t nrTimeIntervals, size_t nrTtotIntervals, bool directedEmissions) {
    auto marginalIndex = new MarginalTripleIndex(nrTimeIntervals, subpopLabels);
    auto lambdaVec = new double[marginalIndex.nrMarginals];
    lambdaVec[] = 1.0;
    return new MSMCmodel(mutationRate, recombinationRate, subpopLabels, lambdaVec, nrTimeIntervals, nrTtotIntervals, directedEmissions);
  }
  
  double emissionProb(string alleles, size_t aij, size_t tTotIndex) const {
    auto triple = marginalIndex.getTripleFromIndex(aij);
    auto id = emissionRate.getEmissionId(alleles, triple.ind1, triple.ind2);
    auto time = timeIntervals.meanTime(triple.time, nrHaplotypes);
    auto tTot = tTotIntervals.meanTime(tTotIndex, 2);
    return emissionRate.emissionProb(id, time, tTot);
  }
  
  double emissionProbHom(size_t time_index, size_t ttotIndex) const {
    auto time = timeIntervals.meanTime(time_index, nrHaplotypes);
    auto tTot = tTotIntervals.meanTime(ttotIndex, 2);
    return emissionRate.emissionProb(0, time, tTot);
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
