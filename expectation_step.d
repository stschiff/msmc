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

class ExpectationStep {
  
  private const MSMCmodel msmc;
  private PropagationCore propagationCore;
  private const string[] fileNames;
  private size_t hmmStrideWidth;
  private double[][] expectationResult;
  private double logLikelihood;
  
  this(in string[] fileNames, MSMCmodel msmc, size_t hmmStrideWidth,
       size_t maxDistance=1000, bool naiveImplementation=false)
  {
    enforce(fileNames.length > 0, "need at least one file");
    enforce(hmmStrideWidth > 0, "need positive hmm stride width");
    this.msmc = msmc;
    this.fileNames = fileNames.dup;
    this.hmmStrideWidth = hmmStrideWidth;
    stderr.writeln("building propagation core");
    if(naiveImplementation)
      propagationCore = new PropagationCoreNaive(msmc, maxDistance);
    else
      propagationCore = new PropagationCoreFast(msmc, maxDistance);
    stderr.writeln("propagation core ready");
    initResults();
  }
  
  private void initResults() {
    expectationResult = new double[][](msmc.nrMarginals, msmc.nrMarginals);
    foreach(ref row; expectationResult)
      row[] = 0.0;
    logLikelihood = 0.0;
  }
  
  void run() {
    foreach(filename; taskPool.parallel(fileNames)) {
      auto result = singleChromosomeExpectation(filename, hmmStrideWidth, propagationCore);
      foreach(au; 0 .. msmc.nrMarginals)
        expectationResult[au][] += result[0][au][];
      logLikelihood += result[1];
    }
  }
  
  private static Tuple!(double[][], double)
  singleChromosomeExpectation(string filename, size_t hmmStrideWidth, in PropagationCore propagationCore)
  {
    auto logTag = format("HMM for file %s: ", filename);
    auto msmc_hmm =  new MSMC_hmm(propagationCore, filename, logTag);

    stderr.writeln(logTag ~ "running forward");
    msmc_hmm.runForward();
    stderr.writeln(logTag ~ "running backward");
    auto exp = msmc_hmm.runBackward(hmmStrideWidth);
    auto logL = msmc_hmm.logLikelihood();
    stderr.writeln(logTag, "likelihood: ", logL);
    msmc_hmm.destroy();
    return tuple(exp, logL);
  }

  private void cleanup() {
    propagationCore.destroy();
    GC.collect();
    GC.minimize();
  }
  
  double[][] getResult() {
    return expectationResult;
  }
  
  double getLogLikelihood() {
    return logLikelihood;
  }

}

unittest {
  writeln("test expectation step");
  auto lambdaVec = new double[12];
  lambdaVec[] = 1.0;
  auto msmc = new MSMCmodel(0.01, 0.001, [0U, 0, 1, 1], lambdaVec, 4, 4);
  auto fileNames = ["model/hmm_testData_tMRCA.txt", "model/hmm_testData_tMRCA.txt", "model/hmm_testData_tMRCA.txt"];
  auto hmmStrideWidth = 100U;
  
  std.parallelism.defaultPoolThreads(1U);
  auto expectationStep = new ExpectationStep(fileNames, msmc, hmmStrideWidth, 100U, true);
  expectationStep.run();
  auto resultSingleThreaded = expectationStep.getResult();
  std.parallelism.defaultPoolThreads(2U);
  expectationStep = new ExpectationStep(fileNames, msmc, hmmStrideWidth, 100U);
  expectationStep.run();
  auto resultMultiThreaded = expectationStep.getResult();
  auto logL = expectationStep.getLogLikelihood();
   
  auto sumSingleThreaded = 0.0;
  auto sumMultiThreaded = 0.0;
  foreach(row; resultSingleThreaded) {
    foreach(val; row) {
      assert(val >= 0.0);
      sumSingleThreaded += val;
    }
  }
  foreach(row; resultMultiThreaded) {
    foreach(val; row) {
      assert(val >= 0.0);
      sumMultiThreaded += val;
    }
  }
  assert(approxEqual(sumSingleThreaded, sumMultiThreaded, 1.0e-8, 0.0), text([sumSingleThreaded, sumMultiThreaded]));
}
