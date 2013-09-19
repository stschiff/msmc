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
import std.math;
import std.string;
import std.conv;
import std.getopt;
import std.parallelism;
import std.algorithm;
import std.array;
import std.json;
import std.file;
import std.typecons;
import std.exception;
import std.c.stdlib;
import core.memory;
import model.msmc_hmm;
import model.msmc_model;
import model.triple_index_marginal;
import model.emission_rate;
import model.transition_rate;
import model.time_intervals;
import model.triple_index_marginal;
import model.coalescence_rate;
import model.rate_integrator;
import model.propagation_core_fastImpl;
import model.data;

void estimateTotalBranchlengths(SegSite_t[] inputData, MSMCmodel params) {

  auto propagationCore = buildPropagationCore(params);
  auto msmc_hmm = buildHMM(inputData, params.nrHaplotypes, propagationCore);
    
  msmc_hmm.runForward();
  
  auto forwardState = propagationCore.newForwardState();
  auto backwardState = propagationCore.newBackwardState();

  foreach_reverse(ref data; inputData) {
    msmc_hmm.getForwardState(forwardState, data.pos);
    msmc_hmm.getBackwardState(backwardState, data.pos);
    double tLeaf = 2.0 * propagationCore.msmc.timeIntervals.meanTime(0, 2);
    auto max = forwardState.vec[0] * backwardState.vec[0];
    foreach(i; 0 .. propagationCore.msmc.nrTimeIntervals) {
      auto p = forwardState.vec[i] * backwardState.vec[i];
      if(p > max) {
        max = p;
        tLeaf = 2.0 * propagationCore.msmc.timeIntervals.meanTime(i, 2);
      }
    }
    data.i_Ttot = params.tTotIntervals.findIntervalForTime(tLeaf);
  }
}
  
private PropagationCoreFast buildPropagationCore(MSMCmodel params) {
  auto lambdaVec = new double[params.nrTtotIntervals];
  lambdaVec[] = 1.0;
  auto expectedTtot =
      params.emissionRate.directedEmissions ? 2.0 : 2.0 * (1.0 + 1.0 / (params.nrHaplotypes - 1.0));
  auto boundaries = TimeIntervals.getQuantileBoundaries(params.nrTtotIntervals, expectedTtot / 2.0);
  auto model = new MSMCmodel(params.mutationRate, params.recombinationRate, [0UL, 0], lambdaVec, boundaries[0 .. $ - 1], 1, params.emissionRate.directedEmissions);

  auto propagationCore = new PropagationCoreFast(model, 1000);
  return propagationCore;
}
  
private MSMC_hmm buildHMM(SegSite_t[] inputData, size_t nrHaplotypes, PropagationCoreFast propagationCore) {
  SegSite_t[] dummyInputData;
  auto alleles = canonicalAlleleOrder(nrHaplotypes);
  foreach(s; inputData) {
    auto dummySite = s.dup;
    if(s.obs.any!"a>1"()) {
      auto count_0 = count(alleles[s.obs[0] - 1], '0');
      auto count_1 = nrHaplotypes - count_0;
      if(propagationCore.msmc.emissionRate.directedEmissions) {
        if(count_1 == 1)
          dummySite.obs = [2];
        else
          dummySite.obs = [1];
      }
      else {
        if(count_0 == 1 || count_1 == 1)
          dummySite.obs = [2];
        else
          dummySite.obs = [1];
      }
    }
    dummyInputData ~= dummySite;
  }
    
  return new MSMC_hmm(propagationCore, dummyInputData);
}

void readTotalBranchlengths(SegSite_t[] inputData, MSMCmodel params, string treeFileName) {
  auto simTreeParser = new SimTreeParser(treeFileName);
  auto allele_order = canonicalAlleleOrder(params.nrHaplotypes);
  foreach(ref segsite; inputData) {
    auto alleles = allele_order[segsite.obs[0] - 1];
    auto tTot = simTreeParser.getTLeafTot(segsite.pos);
    segsite.i_Ttot = params.tTotIntervals.findIntervalForTime(tTot);
  }
}

size_t findDerivedPositition(string alleles) {
  auto pos = 0UL;
  foreach(i, a; alleles) {
    if(a == '1') {
      pos = i;
      break;
    }
  }
  return pos;
}

unittest {
  assert(findDerivedPositition("000100") == 3);
  assert(findDerivedPositition("100000") == 0);
  assert(findDerivedPositition("000001") == 5);
}


class SimTreeParser {
  
  Tuple!(size_t, double)[] data;
  size_t lastIndex;
  
  this(string treeFileName) {
    
    auto treeFile = File(treeFileName, "r");
    auto pos = 0UL;
    foreach(line; treeFile.byLine) {
      auto fields = line.strip().split();
      auto l = fields[0].to!size_t;
      auto str = fields[1];
      auto tLeaflength = getTotLeafLength(str);
      pos += l;
      data ~= tuple(pos, tLeaflength);
    }
  }
  
  double getTLeafTot(size_t pos) {
    auto index = getIndex(pos);
    return data[index][1];
  }
  
  private size_t getIndex(size_t pos) {
    while(data[lastIndex][0] < pos)
      lastIndex += 1;
    while(lastIndex > 0 && data[lastIndex - 1][0] >= pos)
      lastIndex -= 1;
    return lastIndex;
  }
}

double getTotLeafLength(in char[] str) {
  static auto tTotRegex = regex(r"\d+:([\d+\.e-]+)", "g");

  auto matches = match(str, tTotRegex);
  auto times = matches.map!(m => m.captures[1].to!double());
  auto sum = 2.0 * times.reduce!"a+b"();
  return sum;
}

unittest {
  auto tree = "(5:0.1,((2:8.1,1:8.1):0.1,((4:0.1,0:0.1):0.004,3:0.1):0.004):0.01);";
  assert(approxEqual(getTotLeafLength(tree), 2.0 * 16.6, 0.0001, 0.0));
  tree = "((((2:8.3,1:8.3):0.122683,(0:0.11,3:0.11):0.00415405):0.00462688,4:0.12):1.06837,5:1.19);";
  assert(approxEqual(getTotLeafLength(tree), 2.0 * 18.13, 0.0001, 0.0));
}

double[] getLeafLengths(in char[] str) {
  static auto tTotRegex = regex(r"(\d+):([\d+\.e-]+)", "g");

  double[] ret;
  auto matches = match(str, tTotRegex);
  foreach(match; matches) {
    auto i = match.captures[1].to!size_t();
    auto t = 2.0 * match.captures[2].to!double();
    if(i >= ret.length)
      ret.length = i + 1;
    ret[i] = t;
  }
  return ret;
}

unittest {
  auto tree = "((((2:8.3,1:8.3):0.122683,(0:0.11,3:0.11):0.00415405):0.00462688,4:0.12):1.06837,5:1.19);";
  auto leafLengths = getLeafLengths(tree);
  assert(approxEqual(leafLengths[0], 0.22));
  assert(approxEqual(leafLengths[1], 16.6));
  assert(approxEqual(leafLengths[2], 16.6));
  assert(approxEqual(leafLengths[3], 0.22));
  assert(approxEqual(leafLengths[4], 0.24));
  assert(approxEqual(leafLengths[5], 2.38));
  assert(approxEqual(leafLengths.reduce!"a+b"(), getTotLeafLength(tree), 1e-8, 0.0));
}
