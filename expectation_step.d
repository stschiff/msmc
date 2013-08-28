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
 
import std.typecons;
import std.stdio;
import std.string;
import std.algorithm;
import std.exception;
import std.conv;
import std.parallelism;
import std.math;
import core.memory;
import model.propagation_core;
import model.propagation_core_fastImpl;
import model.propagation_core_naiveImpl;
import model.msmc_model;
import model.msmc_hmm;
import model.data;
import logger;

alias Tuple!(double[][], double) ExpectationResult_t;

ExpectationResult_t getExpectation(in SegSite_t[][] inputData, MSMCmodel msmc, size_t hmmStrideWidth,
                                   size_t maxDistance, bool naiveImplementation=false)
{
  PropagationCore propagationCore;
  if(naiveImplementation)
    propagationCore = new PropagationCoreNaive(msmc, maxDistance);
  else
    propagationCore = new PropagationCoreFast(msmc, maxDistance);
  
  auto expectationResult = new double[][](msmc.nrMarginals, msmc.nrMarginals);
  foreach(ref row; expectationResult)
    row[] = 0.0;
  auto logLikelihood = 0.0;
  
  auto cnt = 0;
  foreach(data; taskPool.parallel(inputData)) {
    logInfo(format("\r  * [%s/%s] Expectation Step", ++cnt, inputData.length));
    auto result = singleChromosomeExpectation(data, hmmStrideWidth, propagationCore);
    foreach(au; 0 .. msmc.nrMarginals)
      expectationResult[au][] += result[0][au][];
    logLikelihood += result[1];
  }
  logInfo("\n");
  
  return tuple(expectationResult, logLikelihood);
}

ExpectationResult_t singleChromosomeExpectation(in SegSite_t[] data, size_t hmmStrideWidth,
                                                in PropagationCore propagationCore)
{
  auto msmc_hmm =  new MSMC_hmm(propagationCore, data);

  msmc_hmm.runForward();
  auto exp = msmc_hmm.runBackward(hmmStrideWidth);
  auto logL = msmc_hmm.logLikelihood();
  msmc_hmm.destroy();
  return tuple(exp, logL);
}

unittest {
  writeln("test expectation step");
  auto lambdaVec = new double[12];
  lambdaVec[] = 1.0;
  auto msmc = new MSMCmodel(0.01, 0.001, [0U, 0, 1, 1], lambdaVec, 4, 4);
  auto fileName = "model/hmm_testData.txt";
  auto data = readSegSites(fileName);
  auto hmmStrideWidth = 100UL;
  
  auto allData = [data, data, data];
  
  std.parallelism.defaultPoolThreads(1U);
  auto resultSingleThreaded = getExpectation(allData, msmc, hmmStrideWidth, 100UL, true);
  std.parallelism.defaultPoolThreads(2U);
  auto resultMultiThreaded = getExpectation(allData, msmc, hmmStrideWidth, 100UL);
   
  auto sumSingleThreaded = 0.0;
  auto sumMultiThreaded = 0.0;
  foreach(row; resultSingleThreaded[0]) {
    foreach(val; row) {
      assert(val >= 0.0);
      sumSingleThreaded += val;
    }
  }
  foreach(row; resultMultiThreaded[0]) {
    foreach(val; row) {
      assert(val >= 0.0);
      sumMultiThreaded += val;
    }
  }
  assert(approxEqual(sumSingleThreaded, sumMultiThreaded, 1.0e-8, 0.0), text([sumSingleThreaded, sumMultiThreaded]));
}
