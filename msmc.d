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
import std.file;
import std.typecons;
import std.regex;
import std.exception;
import std.c.stdlib;
import model.data;
import model.msmc_model;
import expectation_step;
import maximization_step;
import logger;
import branchlength;

auto maxIterations = 20UL;
double mutationRate;
double recombinationRate;
size_t[] subpopLabels;
auto timeSegmentPattern = [1UL, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
uint nrThreads;
auto nrTtotSegments = 10UL;
auto verbose = false;
string outFilePrefix;
auto memory = false;
auto naiveImplementation = false;
auto fixedPopSize = false;
auto fixedRecombination = false;
string[] inputFileNames;
size_t hmmStrideWidth = 1000;

auto helpString = "Usage: msmc [options] <datafiles>
  Options:
    -i, --maxIterations=<size_t> : number of EM-iterations [default=20]
    -o, --outFilePrefix=<string> : file prefix to use for all output files
    -m, --mutationRate=<double> : mutation rate, scaled by 2N. In case of more than two haplotypes, this needs to be 
          the same as was used in running \"msmc branchlength\".
    -r, --recombinationRate=<double> : recombination rate, scaled by 2N, to begin with
          [by default set to mutationRate / 4]. Note that recombination rate inference does not behave very well for 
          more than two haplotypes. Using the -R option is recommended for more than 2 haplotypes.
    -t, --nrThreads=<size_t> : nr of threads to use (defaults to nr of CPUs)
    -p, --timeSegmentPattern=<string> : pattern of fixed time segments [default=10*1+15*2]
    -T, --nrTtotSegments=<size_t> : number of discrete values of the total branchlength Ttot [default=10]
    -P, --subpopLabels=<string> comma-separated subpopulation labels (assume one single population by default, with 
          number of haplotypes inferred from first input file). For cross-population analysis with 4 haplotypes, 2 
          coming from each subpopulation, set this to 0,0,1,1
    -R, --fixedRecombination : keep recombination rate fixed [not set by default]
    
    Debug Options:
    -v, --verbose: write out also the transition matrices
    --naiveImplementation: use naive HMM implementation
    --fixedPopSize: learn only the cross-population coalescence rates, keep the population sizes fixed
    --hmmStrideWidth <int> : stride width to traverse the data in the expectation step [=1000]";

void main(string[] args) {
  try {
    parseCommandLine(args);
  }
  catch(Exception e) {
    stderr.writeln("error in parsing command line: ", e.msg);
    exit(0);
  }
  run();
}

void parseCommandLine(string[] args) {
  
  void displayHelpMessageAndExit() {
    stderr.writeln(helpString);
    exit(0);
  }
  void handleTimeSegmentPatternString(string option, string patternString) {
    enforce(match(patternString, r"^\d+\*\d+[\+\d+\*\d+]*"), text("illegal timeSegmentPattern: ", patternString));
    timeSegmentPattern.length = 0;
    foreach(product; std.string.split(patternString, "+")) {
      auto pair = array(map!"to!size_t(a)"(std.string.split(product, "*")));
      foreach(i; 0 .. pair[0]) {
        timeSegmentPattern ~= pair[1];
      }
    }
  }

  void handleSubpopLabelsString(string option, string subpopLabelsString) {
    enforce(match(subpopLabelsString, r"^\d+[,\d+]+"), text("illegal array string: ", subpopLabelsString));
    auto splitted = std.string.split(subpopLabelsString, ",");
    subpopLabels = map!"to!size_t(a)"(splitted).array();
  }
  
  if(args.length == 1) {
    displayHelpMessageAndExit();
  }

  getopt(args,
      std.getopt.config.caseSensitive,
      "maxIterations|i", &maxIterations,
      "mutationRate|m", &mutationRate,
      "recombinationRate|r", &recombinationRate,
      "subpopLabels|P", &handleSubpopLabelsString,
      "timeSegmentPattern|p", &handleTimeSegmentPatternString,
      "nrThreads|t", &nrThreads,
      "nrTtotSegments|T", &nrTtotSegments,
      "verbose|v", &verbose,
      "outFilePrefix|o", &outFilePrefix,
      "help|h", &displayHelpMessageAndExit,
      "naiveImplementation", &naiveImplementation,
      "hmmStrideWidth", &hmmStrideWidth,
      "fixedPopSize", &fixedPopSize,
      "fixedRecombination|R", &fixedRecombination
  );
  if(nrThreads)
    std.parallelism.defaultPoolThreads(nrThreads);
  enforce(!isNaN(mutationRate), "need to set mutation rate");
  if(isNaN(recombinationRate))
    recombinationRate = mutationRate / 4.0;
  enforce(args.length > 1, "need at least one input file");
  enforce(hmmStrideWidth > 0, "hmmStrideWidth must be positive");
  inputFileNames = args[1..$];
  if(subpopLabels.length == 0)
    inferDefaultSubpopLabels();
  
  auto logFileName = outFilePrefix ~ ".log";
  logger.logFile = File(logFileName, "w");
  
  printGlobalParams();
}

void printGlobalParams() {
  logInfo(format("maxIterations:       %s\n", maxIterations));
  logInfo(format("mutationRate:        %s\n", mutationRate));
  logInfo(format("recombinationRate:   %s\n", recombinationRate));
  logInfo(format("subpopLabels:        %s\n", subpopLabels));
  logInfo(format("timeSegmentPattern:  %s\n", timeSegmentPattern));
  logInfo(format("nrThreads:           %s\n", nrThreads));
  logInfo(format("nrTtotSegments:      %s\n", nrTtotSegments));
  logInfo(format("verbose:             %s\n", verbose));
  logInfo(format("outFilePrefix:       %s\n", outFilePrefix));
  logInfo(format("naiveImplementation: %s\n", naiveImplementation));
  logInfo(format("hmmStrideWidth:      %s\n", hmmStrideWidth));
  logInfo(format("fixedPopSize:        %s\n", fixedPopSize));
  logInfo(format("fixedRecombination:  %s\n", fixedRecombination));
}

void inferDefaultSubpopLabels() {
  auto nrHaplotypes = getNrHaplotypesFromFile(inputFileNames[0]);
  foreach(i; 0 .. nrHaplotypes)
    subpopLabels ~= 0;
}

void run() {
  auto nrTimeSegments = reduce!"a+b"(timeSegmentPattern);
  auto params = MSMCmodel.withTrivialLambda(mutationRate, recombinationRate, subpopLabels, nrTimeSegments, 
                                            nrTtotSegments);
  
  auto inputData = inputFileNames.map!(f => readSegSites(f)).array();
  
  
  auto cnt = 0;
  foreach(data; taskPool.parallel(inputData)) {
    logInfo(format("[%s/%s] estimating total branchlengths\r", ++cnt, inputData.length));
    estimateTotalBranchlengths(data, params);
  }
  logInfo("\n");
  
  auto loopFileName = outFilePrefix ~ ".loops.txt";
  
  foreach(iteration; 0 .. maxIterations) {
    logInfo(format("\n[%s/%s] Baumwelch iteration\n", iteration, maxIterations));
    auto expectationResult = getExpectation(inputData, params, hmmStrideWidth, 1000, naiveImplementation);
    auto eMat = expectationResult[0];
    auto logLikelihood = expectationResult[1];
    if(verbose) {
      auto filename = outFilePrefix ~ format(".loop_%s.expectationMatrix.txt", iteration);
      printMatrix(filename, eMat);
    }
    auto newParams = getMaximization(eMat, params, timeSegmentPattern, fixedPopSize, fixedRecombination);
    printLoop(loopFileName, newParams, logLikelihood);
    params = newParams;
  }
  
  auto finalName = outFilePrefix ~ ".final.txt";
  printFinal(finalName, params);
}

void printMatrix(string filename, double[][] eMat) {
  auto f = File(filename, "w");
  foreach(row; eMat) {
    foreach(val; row) {
      f.writef("%s\t", val);
    }
    f.write("\n");
  }
}

void printLoop(string filename, MSMCmodel params, double logLikelihood) {
  auto f = File(filename, "a");
  f.writefln("%s\t%s\t%s", params.recombinationRate, logLikelihood, params.lambdaVec.map!"text(a)".join(",").array());
}

void printFinal(string filename, MSMCmodel params) {
  auto f = File(filename, "w");
  f.write("time_index\ttime_boundary");
  auto nrSubpopPairs = params.nrSubpopulations * (params.nrSubpopulations + 1) / 2;
  foreach(i; 0 .. params.nrSubpopulations) {
    foreach(j; i .. params.nrSubpopulations) {
      f.writef("\tlambda_%s%s", i, j);
    }
  }
  f.write("\n");
  auto lambdaIndex = 0;
  foreach(i; 0 .. params.nrTimeIntervals) {
    f.writef("%s\t%s", i, params.timeIntervals.rightBoundary(i));
    foreach(j; 0 .. nrSubpopPairs) {
      f.writef("\t%s", params.lambdaVec[lambdaIndex++]);
    }
    f.write("\n");
  }
}