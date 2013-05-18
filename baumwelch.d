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
 
import std.json;
import std.exception;
import std.algorithm;
import std.stdio;
import std.file;
import utils;
import model.msmc_model;
import expectation_step;
import maximization_step;

class BaumWelchStepper {
  size_t nrSubThreads;
  bool verbose;
  const size_t[] timeSegmentPattern;
  const string[] inputFileNames;
  size_t hmmStrideWidth;

  MSMCmodel currentIterationParams;
  double logLikelihood;
  ExpectationStep expectationStep;
  
  double[][] expectationResult;
  double[][] modelTransitions;
  
  size_t nrIterationsRun;
  bool naiveImplementation;
  JSONValue report;  
  bool finished;
  bool fixedPopSize;
  bool fixedRecombination;

  this(MSMCmodel initialParams, size_t nrSubThreads, size_t[] timeSegmentPattern, string[] inputFileNames) {
    this.nrSubThreads = nrSubThreads;
    verbose = false;
    this.timeSegmentPattern = timeSegmentPattern;
    this.inputFileNames = inputFileNames;
    naiveImplementation = false;
    
    nrIterationsRun = 0;
    currentIterationParams = initialParams;
    enforce(initialParams.nrTimeIntervals == reduce!"a+b"(timeSegmentPattern.dup));
    logLikelihood = -double.infinity;
    hmmStrideWidth = 1000;
    
    report = JSONValue();
    report.type = JSON_TYPE.ARRAY;
    report.array.length = 0;
    
    finished = false;
  }
  
  void runIteration() {
    stderr.writeln("running Baum-Welch Iteration ", nrIterationsRun + 1);
    stderr.writeln("model parameters:");
    stderr.writeln(currentIterationParams);

    if(nrIterationsRun == 0) {
      expectationStep = new ExpectationStep(nrSubThreads, inputFileNames, currentIterationParams, hmmStrideWidth, 1000, naiveImplementation);
      expectationStep.run();
      logLikelihood = expectationStep.getLogLikelihood();
      stderr.writeln("first logLikelihood: ", logLikelihood);
    }
    
    for(;;) {
      expectationResult = expectationStep.getResult();
      auto maximizationStep = new MaximizationStep(expectationResult, currentIterationParams, timeSegmentPattern);
      maximizationStep.fixedPopSize = fixedPopSize;
      maximizationStep.fixedRecombination = fixedRecombination;
      maximizationStep.run(true);
      auto updatedParams = maximizationStep.getUpdatedParams();
      auto newExpectationStep = new ExpectationStep(nrSubThreads, inputFileNames, updatedParams, hmmStrideWidth, 1000, naiveImplementation);
      newExpectationStep.run();
      // if(newExpectationStep.getLogLikelihood() > logLikelihood) {
        stderr.writeln("found new parameters: ", updatedParams);
        stderr.writefln("previous logLikelihood: %s. new log likelihood: %s", logLikelihood, newExpectationStep.getLogLikelihood());
        currentIterationParams = updatedParams;
        logLikelihood = newExpectationStep.getLogLikelihood();
        modelTransitions = maximizationStep.getModelTransitions();
        expectationStep = newExpectationStep;
        break;
      // }
      // else {
      //   stderr.writefln("previous logLikelihood: %s. new log likelihood: %s", logLikelihood, newExpectationStep.getLogLikelihood());
      //   hmmStrideWidth /= 2;
      //   if(hmmStrideWidth < 50) {
      //     finished = true;
      //     break;
      //   }
      //   stderr.writeln("bad update, trying again with new stride width: ", hmmStrideWidth);
      //   expectationStep = new ExpectationStep(nrSubThreads, inputFileNames, currentIterationParams, hmmStrideWidth, 1000, naiveImplementation);
      //   expectationStep.run();
      // }
    }    
    if(!finished) {
      addToReport();
      nrIterationsRun += 1;
    }
  }

  void addToReport() {
    auto node = JSONValue();
    node.type = JSON_TYPE.OBJECT;
    node.object["updatedParams"] = currentIterationParams.toJSON();
    node.object["hmmStrideWidth"] = makeJSON(hmmStrideWidth);
    node.object["logLikelihood"] = makeJSON(logLikelihood);
    if(verbose) {
      node.object["expectationResult"] = makeJSON(expectationResult);
      node.object["modelTransitions"] = makeJSON(modelTransitions);
    }
    report.array ~= node;
    
  }
  
  JSONValue getReport() {
    return report;
  }
  
}

