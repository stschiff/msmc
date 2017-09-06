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
import core.stdc.stdlib;
import std.range;
import model.data;
import model.msmc_model;
import expectation_step;
import maximization_step;
import logger;
import branchlength;
import model.triple_index_marginal;

auto maxIterations = 20UL;
double mutationRate;
double recombinationRate;
size_t[] subpopLabels;
auto timeSegmentPattern = [1UL, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
uint nrThreads;
size_t nrTtotSegments;
auto verbose = false;
string outFilePrefix;
auto memory = false;
auto naiveImplementation = false;
auto fixedPopSize = false;
auto fixedRecombination = false;
bool directedEmissions = false;
bool skipAmbiguous = false;
string[] inputFileNames, treeFileNames;
SegSite_t[][] inputData;
size_t hmmStrideWidth = 1000;
double[] lambdaVec;
size_t nrTimeSegments;
size_t[] indices;
string logFileName, loopFileName, finalFileName;

immutable versionString = "1.0.1";

auto helpString = format("This is MSMC Version %s. Usage: msmc [options] <datafiles>
  Options:
    -i, --maxIterations=<size_t> : number of EM-iterations [default=20]
    -o, --outFilePrefix=<string> : file prefix to use for all output files
    -r, --rhoOverMu=<double>: ratio of recombination rate over mutation rate. Default=0.25.
    -t, --nrThreads=<size_t> : nr of threads to use (defaults to nr of CPUs)
    -p, --timeSegmentPattern=<string> : pattern of fixed time segments [default=10*1+15*2]
    -P, --subpopLabels=<string> comma-separated subpopulation labels (assume one single population by default, with 
          number of haplotypes inferred from first input file). For cross-population analysis with 4 haplotypes, 2 
          coming from each subpopulation, set this to 0,0,1,1
    -R, --fixedRecombination : keep recombination rate fixed [recommended, but not set by default]
    -I, --indices: indices (comma-separated) of alleles in the data file to run over
    -s, --skipAmbiguous: skip sites with ambiguous phasing. Recommended for gene flow analysis
    -h, --help: show this message", versionString);

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
  
  void handleTreeFileNames(string option, string value) {
    treeFileNames = std.string.split(value, ",").array();
  }
  
  void handleIndices(string option, string value) {
    enforce(match(value, r"^[\d.]+[,[\d.]+]+"), text("illegal indices: ", value));
    indices = std.string.split(value, ",").map!"a.to!size_t()"().array();
  }
  
  if(args.length == 1) {
    displayHelpMessageAndExit();
  }
  
  auto rhoOverMu = 0.25;
  getopt(args,
      std.getopt.config.caseSensitive,
      "maxIterations|i", &maxIterations,
      "mutationRate", &mutationRate,
      "rhoOverMu|r", &rhoOverMu,
      "subpopLabels|P", &handleSubpopLabelsString,
      "timeSegmentPattern|p", &handleTimeSegmentPatternString,
      "nrThreads|t", &nrThreads,
      "nrTtotSegments", &nrTtotSegments,
      "verbose", &verbose,
      "outFilePrefix|o", &outFilePrefix,
      "directedEmissions", &directedEmissions,
      "indices|I", &handleIndices,
      "skipAmbiguous|s", &skipAmbiguous,
      "help|h", &displayHelpMessageAndExit,
      "naiveImplementation", &naiveImplementation,
      "hmmStrideWidth", &hmmStrideWidth,
      "fixedPopSize", &fixedPopSize,
      "fixedRecombination|R", &fixedRecombination,
      "initialLambdaVec", &handleLambdaVecString,
      "treeFileNames", &handleTreeFileNames
  );
  if(nrThreads)
    std.parallelism.defaultPoolThreads(nrThreads);
  enforce(args.length > 1, "need at least one input file");
  enforce(hmmStrideWidth > 0, "hmmStrideWidth must be positive");
  inputFileNames = args[1 .. $];
  if(indices.length == 0)
    inferDefaultIndices();
  if(subpopLabels.length == 0)
    inferDefaultSubpopLabels(indices.length);
  enforce(indices.length == subpopLabels.length, "nr haplotypes in subpopLabels and indices must be equal");
  inputData = readDataFromFiles(inputFileNames, directedEmissions, indices, skipAmbiguous);
  if(isNaN(mutationRate)) {
    stderr.write("estimating scaled mutation rate: ");
    mutationRate = getTheta(inputData, indices.length) / 2.0;
    stderr.writeln(mutationRate);
  }
  recombinationRate = mutationRate * rhoOverMu;
  nrTimeSegments = timeSegmentPattern.reduce!"a+b"();
  if(nrTtotSegments == 0)
    nrTtotSegments = nrTimeSegments;
  auto nrSubpops = MarginalTripleIndex.computeNrSubpops(subpopLabels);
  auto nrMarginals = nrTimeSegments * nrSubpops * (nrSubpops + 1) / 2;
  if(lambdaVec.length > 0) {
    // this is necessary because we read in a scaled lambdaVec.
    lambdaVec[] *= mutationRate;
    enforce(lambdaVec.length == nrMarginals, "initialLambdaVec must have correct length");
  }
  enforce(treeFileNames.length == 0 || treeFileNames.length == inputFileNames.length);
  
  logFileName = outFilePrefix ~ ".log";
  loopFileName = outFilePrefix ~ ".loop.txt";
  finalFileName = outFilePrefix ~ ".final.txt";
  logger.logFile = File(logFileName, "w");
  
  printGlobalParams();
}

void printGlobalParams() {
  logInfo(format("Version:             %s\n", versionString));
  logInfo(format("input files:         %s\n", inputFileNames));
  logInfo(format("maxIterations:       %s\n", maxIterations));
  logInfo(format("mutationRate:        %s\n", mutationRate));
  logInfo(format("recombinationRate:   %s\n", recombinationRate));
  logInfo(format("subpopLabels:        %s\n", subpopLabels));
  logInfo(format("timeSegmentPattern:  %s\n", timeSegmentPattern));
  logInfo(format("nrThreads:           %s\n", nrThreads == 0 ? totalCPUs : nrThreads));
  logInfo(format("nrTtotSegments:      %s\n", nrTtotSegments));
  logInfo(format("verbose:             %s\n", verbose));
  logInfo(format("outFilePrefix:       %s\n", outFilePrefix));
  logInfo(format("naiveImplementation: %s\n", naiveImplementation));
  logInfo(format("hmmStrideWidth:      %s\n", hmmStrideWidth));
  logInfo(format("fixedPopSize:        %s\n", fixedPopSize));
  logInfo(format("fixedRecombination:  %s\n", fixedRecombination));
  logInfo(format("initialLambdaVec:    %s\n", lambdaVec));
  logInfo(format("directedEmissions:   %s\n", directedEmissions));
  logInfo(format("skipAmbiguous:       %s\n", skipAmbiguous));
  logInfo(format("indices:             %s\n", indices));
  logInfo(format("logging information written to %s\n", logFileName));
  logInfo(format("loop information written to %s\n", loopFileName));
  logInfo(format("final results written to %s\n", finalFileName));
  if(verbose)
    logInfo(format("transition matrices written to %s.loop_*.expectationMatrix.txt\n", outFilePrefix));
}

void inferDefaultIndices() {
  auto nrHaplotypes = getNrHaplotypesFromFile(inputFileNames[0]);
  indices = iota(nrHaplotypes).array();
}

void inferDefaultSubpopLabels(size_t nrHaplotypes) {
  subpopLabels = new size_t[nrHaplotypes];
  subpopLabels[] = 0;
}


void run() {
  MSMCmodel params;
  if(lambdaVec.length > 0)
    params = new MSMCmodel(mutationRate, recombinationRate, subpopLabels, lambdaVec, nrTimeSegments, nrTtotSegments, 
                           directedEmissions);
  else
    params = MSMCmodel.withTrivialLambda(mutationRate, recombinationRate, subpopLabels, nrTimeSegments, nrTtotSegments,
                                         directedEmissions);
  
  auto nrFiles = inputData.length;
  if(params.nrHaplotypes > 2) {
    if(treeFileNames.length == 0) {
      auto cnt = 0;
      foreach(i, data; taskPool.parallel(inputData)) {
        logInfo(format("\r[%s/%s] estimating total branchlengths", ++cnt, nrFiles));
        estimateTotalBranchlengths(data, params);
      }
    }
    else {
      auto cnt = 0;
      foreach(data; taskPool.parallel(zip(inputData, treeFileNames))) {
        logInfo(format("\r[%s/%s] estimating total branchlengths", ++cnt, nrFiles));
        readTotalBranchlengths(data[0], params, data[1]);
      }
    }
    logInfo("\n");
  }
  
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

SegSite_t[][] readDataFromFiles(string[] filenames, bool directedEmissions, size_t[] indices, bool skipAmbiguous) {
  SegSite_t[][] ret;
  foreach(filename; filenames) {
    auto data = readSegSites(filename, directedEmissions, indices, skipAmbiguous);
    logInfo(format("read %s SNPs from file %s\n", data.length, filename));
    ret ~= data;
  }
  return ret;
}

void printMatrix(string filename, double[] eVec, double[][] eMat) {
  auto f = File(filename, "w");
  f.writeln(eVec.map!"text(a)"().join("\t"));
  foreach(au; 0 .. eVec.length) {
    f.writeln(eMat[au].map!"text(a)"().join("\t"));
  }
}

void printLoop(string filename, MSMCmodel params, double logLikelihood) {
  auto f = File(filename, "a");
  f.writefln("%s\t%.2f\t%s\t%s", params.recombinationRate, logLikelihood, params.timeIntervals.boundaries.map!(a => text(a * params.mutationRate)).join(",").array(), params.lambdaVec.map!(a => text(a / mutationRate)).join(",").array());
}

void printFinal(string filename, MSMCmodel params) { 
  auto f = File(filename, "w");
  f.write("time_index\tleft_time_boundary\tright_time_boundary");
  auto nrSubpopPairs = params.nrSubpopulations * (params.nrSubpopulations + 1) / 2;
  foreach(i; 0 .. params.nrSubpopulations) {
    foreach(j; i .. params.nrSubpopulations) {
      f.writef("\tlambda_%s%s", i, j);
    }
  }
  f.write("\n");
  auto lambdaIndex = 0;
  foreach(i; 0 .. params.nrTimeIntervals) {
    auto left = params.timeIntervals.leftBoundary(i);
    auto right = params.timeIntervals.rightBoundary(i);
    f.writef("%s\t%s\t%s", i, left * mutationRate, right * mutationRate);
    foreach(j; 0 .. nrSubpopPairs) {
      f.writef("\t%s", params.lambdaVec[lambdaIndex++] / mutationRate);
    }
    f.write("\n");
  }
}
