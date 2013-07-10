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
import std.regex;
import std.exception;
import std.c.stdlib;
import core.memory;
import model.msmc_hmm;
import model.propagation_core_fastImpl;
import model.time_intervals;
import model.triple_index_marginal;
import model.data;
import model.msmc_model;
import powell;
import utils;
import baumwelch;



auto maxIterations = 20UL;
double mutationRate;
double recombinationRate;
size_t[] subpopLabels;
auto timeSegmentPattern = [1UL, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
uint nrThreads;
auto nrTtotSegments = 10UL;
auto verbose = false;
auto outFileName = "/dev/null";
auto memory = false;
auto naiveImplementation = false;
string[] demographyFiles;
auto fixedPopSize = false;
auto fixedRecombination = false;
string[] inputFileNames;
size_t hmmStrideWidth = 1000;

auto helpString = "Usage: msmc inference [options] <datafiles>
  Options:
    -i, --maxIterations=<size_t> : number of EM-iterations [default=20]
    -o, --outFileName=<string> : file to write results to, gets updated every iteration, but includes all previous 
          iterations [default=/dev/null]
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
    --memory: output an estimate for the memory needed and exit
    --naiveImplementation: use naive HMM implementation
    --demographyFiles=<string>: initialize with previous msmc results. This setting must be given once for each 
              subpopulation, as specified in the -P setting. This also overrides the set mutationRate (-m), and the 
              recombination rate (-r) and takes it from the first given demography file
    --fixedPopSize: learn only the cross-population coalescence rates, keep the population sizes fixed
    --hmmStrideWidth <int> : stride width to traverse the data in the expectation step [=1000]";

void inferenceMain(string[] args) {
  try {
    parseCommandLine(args);
  }
  catch(Exception e) {
    stderr.writeln("error in parsing command line: ", e.msg);
    exit(0);
  }
  if(subpopLabels.length == 0)
    inferDefaultSubpopLabels();
  
  if(memory)
    runMemory();
  else
    runBaumWelch();
}

void parseCommandLine(string[] args) {
  
  void displayHelpMessageAndExit() {
    stderr.writeln(helpString);
    exit(0);
  }
  void handleTimeSegmentPatternString(string option, string timeSegmentPatternString) {
    timeSegmentPattern = parseTimeSegmentPattern(timeSegmentPatternString);
  }

  void handleSubpopLabelsString(string option, string subpopLabelsString) {
    subpopLabels = parseCommaSeparatedArray(subpopLabelsString);
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
      "outFileName|o", &outFileName,
      "memory", &memory,
      "help|h", &displayHelpMessageAndExit,
      "naiveImplementation", &naiveImplementation,
      "demographyFiles", &demographyFiles,
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
}

void inferDefaultSubpopLabels() {
  auto nrHaplotypes = getNrHaplotypesFromFile(inputFileNames[0]);
  stderr.writeln("found ", nrHaplotypes, " haplotypes in file");
  foreach(i; 0 .. nrHaplotypes)
    subpopLabels ~= 0;
}

void runMemory() {
  
  auto nrTimeSegments = reduce!"a+b"(timeSegmentPattern);
  auto msmc = MSMCmodel.withTrivialLambda(mutationRate, recombinationRate, subpopLabels, nrTimeSegments, nrTtotSegments);
  
  auto maxDistance = 1000;
  int bytesPerDouble = 8;
    
  // plus one for missing data
  auto memoryForDiags = msmc.nrMarginals * (nrTtotSegments + 1) * maxDistance * bytesPerDouble;
  auto memoryForOffDiags = msmc.nrMarginals ^^ 2 * (nrTtotSegments + 1) * maxDistance * bytesPerDouble;
  // forward and backward
  auto memoryForPropagators = memoryForDiags + 2 * memoryForOffDiags;
    
  auto nrFiles = std.parallelism.defaultPoolThreads;
  if(nrFiles > inputFileNames.length)
    nrFiles = cast(uint)inputFileNames.length;
    
  auto ttotDummy = TimeIntervals.standardTotalBranchlengthIntervals(1, 2);
  int totalNrSegSites = 0;
  writefln("reading %s files", nrFiles);
  foreach(i; 0 .. nrFiles) {
    auto segsites = readSegSites(inputFileNames[i], cast(size_t)subpopLabels.length, ttotDummy);
    auto chopped = chop_segsites(segsites, maxDistance);
    totalNrSegSites += cast(int)chopped.length;
  }
  
  writeln("nrStates: ", msmc.nrStates, ", nrMarginals: ", msmc.nrMarginals);
  // scaling vector, forward states, marginal states plus 2 pointers
  auto memoryForStates = cast(long)bytesPerDouble * cast(long)totalNrSegSites * cast(long)(msmc.nrStates + msmc.nrMarginals + 1 + 2);
  // 2 size_ts plus obs-vec : 2 pointers plus 64 bytes initially reserved = 80 bytes
  double memoryForSegSites = totalNrSegSites * 80;
  auto memoryForHmm = memoryForStates + memoryForSegSites;
    
  auto gb = 1073741824.0;
  writeln("nr of subpopulations: ", msmc.nrSubpopulations);
  writefln("memory needed for precomputed propagation matrices: %s GB", memoryForPropagators / gb);
  writeln("total number of positions to store state probabilities: ", totalNrSegSites);
  writefln("total memory needed for state probabilities: %s GB", memoryForHmm / gb);
  writefln("TOTAL: %s GB", (memoryForPropagators + memoryForHmm) / gb);
    
  exit(0);
}

void runBaumWelch() {
  auto nrTimeSegments = reduce!"a+b"(timeSegmentPattern);
  auto msmc = MSMCmodel.withTrivialLambda(mutationRate, recombinationRate, subpopLabels, nrTimeSegments, nrTtotSegments);
  if(demographyFiles)
    msmc = MSMCmodel.overrideDemographies(msmc, demographyFiles);
  
  auto baumWelchStepper = new BaumWelchStepper(msmc, timeSegmentPattern, inputFileNames);
  baumWelchStepper.verbose = verbose;
  baumWelchStepper.naiveImplementation = naiveImplementation;
  baumWelchStepper.fixedPopSize = fixedPopSize;
  baumWelchStepper.fixedRecombination = fixedRecombination;
  baumWelchStepper.hmmStrideWidth = hmmStrideWidth;
  
  while(baumWelchStepper.nrIterationsRun < maxIterations) {
    baumWelchStepper.runIteration();
    if(baumWelchStepper.finished)
      break;
    auto resultReport = baumWelchStepper.getReport();

    JSONValue json;
    json.type = JSON_TYPE.OBJECT;
    json.object["nrIterationsRun"] = makeJSON(baumWelchStepper.nrIterationsRun);
    json.object["initialParams"] = msmc.toJSON();
    json.object["timeSegmentPattern"] = makeJSON(timeSegmentPattern);
    json.object["inputFileNames"] = makeJSON(inputFileNames);
    json.object["results"] = resultReport;
    auto outFile = File(outFileName, "w");
    outFile.write(toJSON(&json));
    outFile.close();
  }
}

