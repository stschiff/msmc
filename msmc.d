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
auto nrTtotInternal = 40UL;
auto verbose = false;
string outFilePrefix;
auto memory = false;
auto naiveImplementation = false;
auto fixedPopSize = false;
auto fixedRecombination = false;
string[] inputFileNames;
size_t hmmStrideWidth = 1000;
double[] lambdaVec;
size_t nrTimeSegments;
string logFileName, loopFileName, finalFileName;


auto helpString = "Usage: msmc [options] <datafiles>
  Options:
    -i, --maxIterations=<size_t> : number of EM-iterations [default=20]

    -o, --outFilePrefix=<string> : file prefix to use for all output files

    -m, --mutationRate=<double> : mutation rate, scaled by 2N. In case of more than two haplotypes, this needs to be 
          the same as was used in running \"msmc branchlength\".

    -r, --recombinationRate=<double> : recombination rate, scaled by 2N, to begin with
          [by default set to mutationRate / 4]. Recombination rate inference does not work very well for 
          more than two haplotypes. Using the -R option is recommended for more than 2 haplotypes.

    -t, --nrThreads=<size_t> : nr of threads to use (defaults to nr of CPUs)

    -p, --timeSegmentPattern=<string> : pattern of fixed time segments [default=10*1+15*2]

    -T, --nrTtotSegments=<size_t> : number of discrete values of the total branchlength Ttot [default=10]

    -O, --nrTtotInternal=<size_t> : number of states used for the Ttot-HMM [default=40]

    -P, --subpopLabels=<string> comma-separated subpopulation labels (assume one single population by default, with 
          number of haplotypes inferred from first input file). For cross-population analysis with 4 haplotypes, 2 
          coming from each subpopulation, set this to 0,0,1,1

    -R, --fixedRecombination : keep recombination rate fixed [recommended, but not set by default]

    -v, --verbose: write out the expected number of transition matrices (into a separate file)

    --naiveImplementation: use naive HMM implementation [for debugging only]

    --fixedPopSize: learn only the cross-population coalescence rates, keep the population sizes fixed [not recommended]

    --hmmStrideWidth <int> : stride width to traverse the data in the expectation step [default=1000]

    --initialLambdaVec <str> : comma-separated string of lambda-values to start with. This can be used to
      continue a previous run by copying the values in the last row and the third column of the corresponding
      *.loop file";

void main(string[] args) {
  try {
    parseCommandLine(args);
  }
  catch(Exception e) {
    stderr.writeln("error in parsing command line: ", e);
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
  
  void handleLambdaVecString(string option, string lambdaString) {
    enforce(match(lambdaString, r"^[\d.]+[,[\d.]+]+"), text("illegal array string: ", lambdaString));
    lambdaVec = std.string.split(lambdaString, ",").map!"to!double(a)"().array();
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
      "nrTtotInternal|O", &nrTtotInternal,
      "verbose|v", &verbose,
      "outFilePrefix|o", &outFilePrefix,
      "help|h", &displayHelpMessageAndExit,
      "naiveImplementation", &naiveImplementation,
      "hmmStrideWidth", &hmmStrideWidth,
      "fixedPopSize", &fixedPopSize,
      "fixedRecombination|R", &fixedRecombination,
      "initialLambdaVec", &handleLambdaVecString
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
  nrTimeSegments = reduce!"a+b"(timeSegmentPattern);
  if(lambdaVec.length == 0) {
    lambdaVec = new double[nrTimeSegments];
    lambdaVec[] = 1.0;
  }
  enforce(lambdaVec.length == nrTimeSegments, "initialLambdaVec must have correct length");
  
  logFileName = outFilePrefix ~ ".log";
  loopFileName = outFilePrefix ~ ".loop.txt";
  finalFileName = outFilePrefix ~ ".final.txt";
  logger.logFile = File(logFileName, "w");
  
  printGlobalParams();
}

void printGlobalParams() {
  logInfo(format("maxIterations:       %s\n", maxIterations));
  logInfo(format("mutationRate:        %s\n", mutationRate));
  logInfo(format("recombinationRate:   %s\n", recombinationRate));
  logInfo(format("subpopLabels:        %s\n", subpopLabels));
  logInfo(format("timeSegmentPattern:  %s\n", timeSegmentPattern));
  logInfo(format("nrThreads:           %s\n", nrThreads == 0 ? totalCPUs : nrThreads));
  logInfo(format("nrTtotSegments:      %s\n", nrTtotSegments));
  logInfo(format("nrTtotInternal:      %s\n", nrTtotInternal));
  logInfo(format("verbose:             %s\n", verbose));
  logInfo(format("outFilePrefix:       %s\n", outFilePrefix));
  logInfo(format("naiveImplementation: %s\n", naiveImplementation));
  logInfo(format("hmmStrideWidth:      %s\n", hmmStrideWidth));
  logInfo(format("fixedPopSize:        %s\n", fixedPopSize));
  logInfo(format("fixedRecombination:  %s\n", fixedRecombination));
  logInfo(format("initialLambdaVec:    %s\n", lambdaVec));
  logInfo(format("logging information written to %s\n", logFileName));
  logInfo(format("loop information written to %s\n", loopFileName));
  logInfo(format("final results written to %s\n", finalFileName));
  if(verbose)
    logInfo(format("transition matrices written to %s.loop_*.expectationMatrix.txt\n", outFilePrefix));
}

void inferDefaultSubpopLabels() {
  auto nrHaplotypes = getNrHaplotypesFromFile(inputFileNames[0]);
  foreach(i; 0 .. nrHaplotypes)
    subpopLabels ~= 0;
}

void run() {
  auto params = new MSMCmodel(mutationRate, recombinationRate, subpopLabels, lambdaVec, nrTimeSegments, nrTtotSegments);
  
  auto inputData = readDataFromFiles(inputFileNames);
  
  auto cnt = 0;
  foreach(data; taskPool.parallel(inputData)) {
    logInfo(format("\r[%s/%s] estimating total branchlengths", ++cnt, inputData.length));
    estimateTotalBranchlengths(data, params, nrTtotInternal);
  }
  logInfo("\n");
  
  auto f = File(loopFileName, "w");
  f.close();
  
  foreach(iteration; 0 .. maxIterations) {
    logInfo(format("[%s/%s] Baumwelch iteration\n", iteration + 1, maxIterations));
    auto expectationResult = getExpectation(inputData, params, hmmStrideWidth, 1000, naiveImplementation);
    auto eVec = expectationResult[0];
    auto eMat = expectationResult[1];
    auto logLikelihood = expectationResult[2];
    printLoop(loopFileName, params, logLikelihood);
    if(verbose) {
      auto filename = outFilePrefix ~ format(".loop_%s.expectationMatrix.txt", iteration);
      printMatrix(filename, eVec, eMat);
    }
    auto newParams = getMaximization(eVec, eMat, params, timeSegmentPattern, fixedPopSize, fixedRecombination);
    params = newParams;
  }
  
  printFinal(finalFileName, params);
}

SegSite_t[][] readDataFromFiles(string[] filenames) {
  SegSite_t[][] ret;
  foreach(filename; filenames) {
    auto data = readSegSites(filename);
    logInfo(format("read %s SNPs from file %s\n", data.length, filename));
    ret ~= data;
  }
  return ret;
}

void printMatrix(string filename, double[] eVec, double[][] eMat) {
  auto f = File(filename, "w");
  foreach(au; 0 .. eVec.length) {
    foreach(bv; 0 .. eVec.length) {
      auto val = eMat[au][bv];
      if(au == bv)
        val += eVec[au];
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