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
 
module model.msmc_hmm;
import std.stdio;
import std.math;
import std.string;
import std.conv;
import std.algorithm;
import std.concurrency;
import std.typecons;
import std.random;
import std.exception;
import model.data;
import model.gsl_matrix_vector;
import model.propagation_core;
import model.stateVec;
import model.stateVecAllocator;

SegSite_t[] chop_segsites(in SegSite_t[] segsites, size_t maxDistance) {
  SegSite_t[] ret;
  size_t lastPos = 0;
  foreach(segsite; segsites) {
    while(segsite.pos - lastPos > maxDistance) {
      ret ~= new SegSite_t(lastPos + maxDistance, min(segsite.obs[0], 1), segsite.i_Ttot);
      lastPos += maxDistance;
    }
    ret ~= segsite.dup;
    lastPos = segsite.pos;
  }
  return ret;
}

unittest {
  writeln("testing chop_sites");
  auto data = [
    new SegSite_t(400, [3], 0),
    new SegSite_t(3600, [0], 0), // missing data
    new SegSite_t(5000, [2], 0)
  ];
  auto ret = chop_segsites(data, 1000);
  assert(ret[0].pos == 400);
  assert(ret[0].obs == [3]);
  assert(ret[1].pos == 1400);
  assert(ret[1].obs == [0]);
  assert(ret[3].pos == 3400);
  assert(ret[3].obs == [0]);
  assert(ret[5].pos == 4600);
  assert(ret[5].obs == [1]); // homozygous
  assert(ret[6].pos == 5000);
  assert(ret[6].obs == [2]);
}

class MSMC_hmm {

  size_t L;
  size_t maxDistance;
  size_t indexCache;
  size_t hmmStrideWidth;
  size_t currentBackwardIndex;
  State_t currentBackwardState, nextBackwardState;
  
  const SegSite_t[] segsites;
  double[] scalingFactors;
  const PropagationCore propagationCore;
  StateVecAllocator stateVecAllocator;
  State_t[] forwardStates;
  State_t expectationForwardDummy, expectationBackwardDummy, getBackwardStateDummy;
  State_t runBackwardDummy;
  bool have_run_forward;
  
  this(in PropagationCore propagationCore, in SegSite_t[] segsites) {
    this.propagationCore = propagationCore;
    this.maxDistance = propagationCore.maxDistance;
    this.segsites = chop_segsites(segsites, maxDistance);
    this.L = this.segsites.length;
    this.hmmStrideWidth = hmmStrideWidth;
    
    scalingFactors = new double[L];
    scalingFactors[] = 0.0;
    
    auto stateSize = propagationCore.forwardStateSize;
    stateVecAllocator = new StateVecAllocator(L * stateSize);
    forwardStates = new State_t[L];
    foreach(i; 0 .. L)
      forwardStates[i] = propagationCore.newForwardState(stateVecAllocator);

    expectationForwardDummy = propagationCore.newForwardState();
    expectationBackwardDummy = propagationCore.newBackwardState();
    getBackwardStateDummy = propagationCore.newBackwardState();
    runBackwardDummy = propagationCore.newBackwardState();
    currentBackwardState = propagationCore.newBackwardState();
    nextBackwardState = propagationCore.newBackwardState();
    currentBackwardIndex = L - 1;
  }

  double logLikelihood() const {
    return scalingFactors.map!log().reduce!"a+b"();
  }
  
  void runForward()
  out {
    foreach(i; 0 .. L)
      assert(scalingFactors[i] > 0);
  }
  body {
    enforce(!have_run_forward);
    propagationCore.initialState(forwardStates[0]);
    scalingFactors[0] = forwardStates[0].norm;
    forwardStates[0].scale(1.0 / scalingFactors[0]);
    

    auto forwardDummyVec = propagationCore.newForwardState();
    foreach(index; 1 .. L) {
      
      if(segsites[index].pos == segsites[index - 1].pos + 1) {
        propagationCore.propagateSingleForward(forwardStates[index - 1],
            forwardStates[index], segsites[index - 1], segsites[index]);
      }
      else {
        auto dummy_site = getSegSite(segsites[index].pos - 1);
        propagationCore.propagateMultiForward(forwardStates[index - 1], 
            forwardDummyVec, segsites[index - 1], dummy_site);
        propagationCore.propagateSingleForward(forwardDummyVec, forwardStates[index],
            dummy_site, segsites[index]);
      }
      scalingFactors[index] = forwardStates[index].norm;
      assert(scalingFactors[index] > 0.0, text(scalingFactors));
      forwardStates[index].scale(1.0 / scalingFactors[index]);
    }
    have_run_forward = true;
  }

  double[][] runBackward(size_t hmmStrideWidth=1000) {
    enforce(have_run_forward);

    auto nrMarginals = propagationCore.getMSMC.nrMarginals;

    auto forwardBackwardResult = new double[][](nrMarginals, nrMarginals);
    foreach(i; 0 .. nrMarginals)
      forwardBackwardResult[i][] = 0.0;
    
    currentBackwardIndex = L - 1;
    auto expec = new double[][](nrMarginals, nrMarginals);
    for(size_t pos = segsites[$ - 1].pos; pos > segsites[0].pos && pos <= segsites[$ - 1].pos; pos -= hmmStrideWidth) {
      getSingleExpectation(pos, expec);
      foreach(i; 0 .. nrMarginals)
        forwardBackwardResult[i][] += expec[i][];
    }

    return forwardBackwardResult;
  }
  
  private void getSingleExpectation(size_t pos, double[][] ret)
  in {
    assert(pos > segsites[0].pos, text(pos, " ", segsites[0].pos));
    assert(pos <= segsites[$ - 1].pos, text([pos, segsites[0].pos]));
    assert(have_run_forward);
  }
  out {
    auto sum = 0.0;
    foreach(row; ret)
      sum += reduce!"a+b"(row);
    assert(approxEqual(sum, 1.0, 1.0e-8, 0.0), text(sum));
  }
  body {    
    getForwardState(expectationForwardDummy, pos - 1);
    getBackwardState(expectationBackwardDummy, pos);
    auto site = getSegSite(pos);
    
    propagationCore.getTransitionExpectation(expectationForwardDummy, expectationBackwardDummy, site, ret);
  } 
  
  void getForwardState(State_t s, size_t pos)
  in {
    assert(pos >= segsites[0].pos);
    assert(pos <= segsites[$ - 1].pos);
  }
  body {
    auto index = getRightIndexAtPos(pos);
    if(pos == segsites[index].pos) {
      forwardStates[index].copy_into(s);
    }
    else {
      auto site = getSegSite(pos);
      propagationCore.propagateMultiForward(forwardStates[index - 1], s, segsites[index - 1], site);
    }
  }

  void getBackwardState(State_t s, size_t pos)
  in {
    assert(pos >= segsites[0].pos);
    assert(pos <= segsites[$ - 1].pos);
  }
  body {
    auto index = getRightIndexAtPos(pos);
    auto site = getSegSite(pos);
    if(pos == segsites[index].pos) {
      getBackwardStateAtIndex(index).copy_into(s);
    }
    else {
      if(pos == segsites[index].pos - 1) {
        propagationCore.propagateSingleBackward(getBackwardStateAtIndex(index), s, segsites[index], site);
      }
      else {
        auto dummy_site = getSegSite(segsites[index].pos - 1);
        propagationCore.propagateSingleBackward(getBackwardStateAtIndex(index), getBackwardStateDummy, segsites[index], dummy_site);
        propagationCore.propagateMultiBackward(getBackwardStateDummy, s, dummy_site, site);
      }
    }
  }
  
  private SegSite_t getSegSite(size_t pos) {
    auto index = getRightIndexAtPos(pos);
    if(segsites[index].pos == pos)
      return segsites[index].dup;
    else
      return new SegSite_t(pos, min(segsites[index].obs[0], 1), segsites[index].i_Ttot);
  }
  
  private size_t getRightIndexAtPos(size_t pos)
  in {
    assert(pos <= segsites[L - 1].pos);
  }
  out(result) {
    assert(segsites[result].pos >= pos);
    if(result > 0) {
      assert(segsites[result - 1].pos < pos);
    }
  }
  body {
    while(segsites[indexCache].pos < pos) {
      ++indexCache;
    }
    if(indexCache > 0) {
      while(segsites[indexCache - 1].pos >= pos) {
        if(--indexCache == 0)
            break;
      }
    }
    return indexCache;
  }

  private State_t getBackwardStateAtIndex(size_t index)
  in {
    assert(have_run_forward);
    assert(index < L);
  }
  body {
    if(index == L - 1) {
      assert(scalingFactors[L - 1] > 0, text(scalingFactors[L - 1]));
      propagationCore.setState(currentBackwardState, 1.0 / scalingFactors[L - 1], segsites[L - 1]);
      currentBackwardIndex = L - 1;
    }
    else {
      assert(index <= currentBackwardIndex, text([index, L]));
      while(index < currentBackwardIndex) {
        if(segsites[currentBackwardIndex].pos == segsites[currentBackwardIndex - 1].pos + 1) {
          propagationCore.propagateSingleBackward(currentBackwardState, nextBackwardState,
              segsites[currentBackwardIndex], segsites[currentBackwardIndex - 1]);
        }
        else {
          auto dummy_site = getSegSite(segsites[currentBackwardIndex].pos - 1);
          propagationCore.propagateSingleBackward(currentBackwardState, runBackwardDummy,
              segsites[currentBackwardIndex], dummy_site);
          propagationCore.propagateMultiBackward(runBackwardDummy, nextBackwardState,
              dummy_site, segsites[currentBackwardIndex - 1]);
        }
        nextBackwardState.copy_into(currentBackwardState);
        currentBackwardState.scale(1.0 / scalingFactors[currentBackwardIndex - 1]);
        currentBackwardIndex -= 1;
      }
    }
    return currentBackwardState;
  }

}

unittest {
  writeln("testing MSMC_hmm");
  import model.propagation_core_naiveImpl;
  import model.propagation_core_fastImpl;
  import model.msmc_model;
  
  auto lambdaVec = new double[30];
  lambdaVec[] = 1.0;
  auto params = new MSMCmodel(0.01, 0.001, [0U, 0, 1, 1], lambdaVec, 10, 4);
  auto lvl = 1.0e-8;
  
  auto propagationCoreNaive = new PropagationCoreNaive(params, 100);
  auto propagationCoreFast = new PropagationCoreFast(params, 100);

  auto nrS = propagationCoreFast.getMSMC.nrStates;
  
  auto data = readSegSites("model/hmm_testData.txt");
  
  auto msmc_hmm_fast = new MSMC_hmm(propagationCoreFast, data);
  auto msmc_hmm_naive = new MSMC_hmm(propagationCoreNaive, data);
  msmc_hmm_fast.runForward();
  msmc_hmm_naive.runForward();
  
  
  auto L = msmc_hmm_naive.L;
  
  foreach(pos; 0 .. msmc_hmm_fast.L) {
    assert(approxEqual(msmc_hmm_fast.scalingFactors[pos],
                       msmc_hmm_naive.scalingFactors[pos], lvl, 0.0),
           text([pos, msmc_hmm_fast.scalingFactors[pos], 
                msmc_hmm_naive.scalingFactors[pos]]));
  }
  
  for(auto pos = L - 1; pos >= 0 && pos < L; --pos) {
    foreach(aij; 0 .. nrS) {
      assert(
          approxEqual(
              msmc_hmm_naive.forwardStates[pos].vec[aij],
              msmc_hmm_fast.forwardStates[pos].vec[aij],
              lvl, 0.0
          ),
          text(
              [msmc_hmm_naive.forwardStates[pos].vec[aij],
              msmc_hmm_fast.forwardStates[pos].vec[aij]]
          )
      );
      assert(
          approxEqual(
              msmc_hmm_naive.getBackwardStateAtIndex(pos).vec[aij],
              msmc_hmm_fast.getBackwardStateAtIndex(pos).vec[aij],
              lvl, 0.0
          ),
          text(
              [pos, msmc_hmm_naive.getBackwardStateAtIndex(pos).vec[aij],
              msmc_hmm_fast.getBackwardStateAtIndex(pos).vec[aij]]
          )
      );
    }
  }
  for(auto pos = L - 1; pos >= 0 && pos < L; --pos) {
    auto sum_f = 0.0;
    auto sum_n = 0.0;
    foreach(aij; 0 .. propagationCoreFast.getMSMC.nrStates) {
      sum_f += msmc_hmm_fast.forwardStates[pos].vec[aij] *
        msmc_hmm_fast.getBackwardStateAtIndex(pos).vec[aij] * 
        msmc_hmm_fast.scalingFactors[pos];
      sum_n += msmc_hmm_naive.forwardStates[pos].vec[aij] *
        msmc_hmm_naive.getBackwardStateAtIndex(pos).vec[aij] * 
        msmc_hmm_naive.scalingFactors[pos];
    }
    
    assert(approxEqual(sum_f, 1.0, lvl, 0.0), text(sum_f));
    assert(approxEqual(sum_n, 1.0, lvl, 0.0), text(sum_n));
  }
  auto expec = msmc_hmm_fast.runBackward();
  auto expec_n = msmc_hmm_naive.runBackward();
  foreach(alpha; 0 .. params.nrTimeIntervals) {
    foreach(beta; 0 .. params.nrTimeIntervals) {
      assert(
          approxEqual(expec[alpha][beta], expec_n[alpha][beta], lvl, 0.0),
          text([expec[alpha][beta], expec_n[alpha][beta]])
      );
    }
  }
}


