import std.stdio;
import std.math;
import std.string;
import std.conv;
import std.getopt;
import std.concurrency;
import std.algorithm;
import std.array;
import std.json;
import std.file;
import std.typecons;
import std.regex;
import std.exception;
import std.c.stdlib;
import core.memory;
import msmc_hmm;
import utils;
import msmc_utils;
import msmc_model;
import triple_index_marginal;
import expectation_step;
import time_intervals;
import transition_rate;
import data;

class ExpectationApplication {
  
  double mutationRate;
  size_t nrSubThreads;
  size_t nrTtotSegments;
  size_t[] subpopLabels;
  string[] inputFileNames;
  size_t nrTimeSegments;
  string[] demographyFiles;
  double recombinationRate;
  
  string subpopLabelsString;
  
  static void runWithArgs(string[] args) {
    auto app = new ExpectationApplication(args);
    app.run();
  }
  
  this(string[] args) {
    nrTimeSegments = 40;
    nrSubThreads = 1;
    nrTtotSegments = 10;
    parseCommandLine(args);
  }

  void parseCommandLine(string[] args) {
    if(args.length == 1)
      displayHelpMessageAndExit();
    try {
      readArguments(args);
      if(subpopLabelsString)
        subpopLabels = parseCommaSeparatedArray(subpopLabelsString);
      else {
        auto nrHaplotypes = getNrHaplotypesFromFile(inputFileNames[0]);
        stderr.writeln("found ", nrHaplotypes, " haplotypes in file");
        foreach(i; 0 .. nrHaplotypes)
          subpopLabels ~= 0;
      }
      validate();
    }
    catch(Exception e) {
      displayHelpMessageAndExit(e);
    }
  }
  
  void readArguments(string[] args) {
    getopt(args,
        std.getopt.config.caseSensitive,
        "mutmationRate|m", &mutationRate,
        "recombinationRate|r", &recombinationRate,
        "subpopLabels|P", &subpopLabelsString,
        "nrSubThreads|t", &nrSubThreads,
        "nrTtotSegments|T", &nrTtotSegments,
        "nrTimeSegments|n", &nrTimeSegments,
        "demographyFiles", &demographyFiles
    );
    inputFileNames = args[1..$];
  }

  private void validate() {
    enforce(inputFileNames.length > 0, "need at least one input file");
  }
  
  static void displayHelpMessageAndExit() {
    stderr.writeln("Usage: msmc expectation [options] <datafiles>
    Options:

      -n, --nrTimeSegments=<int> : nr of time intervals [=40]

      -m, --mutationRate=<double> : scaled mutation rate

      -r, --recombinationRate=<double> : scaled recombination rate

      -t, --nrSubThreads=<size_t> : nr of threads [=1]

      -T, --nrTtotSegments=<size_t> : number of discrete values of Ttot [=10]

      -P, --subpopLabels=<string> comma-separated subpopulation labels [optional]

      --demographyFiles=<string>: begin iterations with a given demography and recombination rate from a file. This 
                settings may be needed multiple times, once for each subpopulations, as specified in the -P setting");
    exit(0);
  }
  
  static void displayHelpMessageAndExit(Exception e) {
    writeln(e.msg);
    displayHelpMessageAndExit();
  }
  
  void run() {
    auto msmc = MSMCmodel.withTrivialLambda(mutationRate, recombinationRate, subpopLabels, nrTimeSegments, nrTtotSegments);
    if(demographyFiles)
      msmc = MSMCmodel.overrideDemographies(msmc, demographyFiles);
    
    ExpectationStep expectationStep;
    try {
      expectationStep = new ExpectationStep(nrSubThreads, inputFileNames, msmc, 1000);
    }
    catch(Exception e) {
      displayHelpMessageAndExit(e);
    }
    
    expectationStep.run();
    auto result = expectationStep.getResult();
    writeMatrix(result);    
    auto modelTransitions = computeModelTransitions(msmc);
    writeMatrix(modelTransitions);
    
  }
  
  double[][] computeModelTransitions(MSMCmodel msmc) {
    auto mat = new double[][](msmc.nrMarginals, msmc.nrMarginals);
    foreach(au; 0 .. msmc.nrMarginals) {
      foreach(bv; 0 .. msmc.nrMarginals) {
        mat[au][bv] = msmc.transitionRate.transitionProbabilityMarginal(au, bv);
      }
    }
    return mat;
  }
  
  void writeMatrix(double[][] matrix) {
    foreach(row; matrix) {
      foreach(val; row)
        write(val, " ");
      writeln("");
    }
  }
}  
