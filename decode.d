import std.stdio;
import std.getopt;
import std.exception;
import std.c.stdlib;
import std.algorithm;
import msmc_hmm;
import data;
import time_intervals;
import msmc_model;
import propagation_core_fastImpl;

double mutationRate, recombinationRate;
size_t nrTimeSegments=40, nrTtotSegments=40, stride=1000;
bool tTot = false;
string inputFileName;
size_t nrHaplotypes;

void decodeMain(string[] args) {
  try {
    parseCommandlineArgs(args);
  }
  catch (Exception e) {
    stderr.writeln(e.msg);
    displayHelpMessage();
    exit(0);
  }
  run();
}

void parseCommandlineArgs(string[] args) {
  getopt(args,
         std.getopt.config.caseSensitive,
         "mutationRate|m", &mutationRate,
         "recombinationRate|r", &recombinationRate,
         "nrTimeSegments|t", &nrTimeSegments,
         "nrTtotSegments|T", &nrTtotSegments,
         "tTot", &tTot,
         "stride|s", &stride);
  enforce(args.length == 2, "need exactly one input file");
  inputFileName = args[1];
  nrHaplotypes = getNrHaplotypesFromFile(inputFileName);
  enforce(mutationRate > 0, "need positive mutationRate");
  enforce(recombinationRate > 0, "need positive recombinationRate");
}

void displayHelpMessage() {
  stderr.writeln("Usage: msmc decode [options] <inputFile>
Options:
-m, --mutationRate <double>
-r, --recombinationRate <double>
-t, --nrTimeSegments <int>
-T, --nrTtotSegments <int>
--tTot : output the total branch length rather than the first coalescent time
-s, --stride <int>");
}

void run() {
  auto hmm = tTot ? makeTtotHmm() : makeStandardHmm();
  decodeWithHmm(hmm);
}

MSMC_hmm makeTtotHmm() {
  auto lambdaVec = new double[nrTtotSegments];
  lambdaVec[] = 1.0;
  auto expectedTtot = 2.0 * TimeIntervals.computeWattersonFactor(nrHaplotypes);
  auto boundaries = TimeIntervals.getQuantileBoundaries(nrTtotSegments, expectedTtot / 2.0);
  auto model = new MSMCmodel(mutationRate, recombinationRate, [0UL, 0], lambdaVec, boundaries[0 .. $ - 1], 1);

  stderr.writeln("generating propagation core");
  auto propagationCore = new PropagationCoreFast(model, 1000);

  stderr.writeln("loading file ", inputFileName);
  auto segsites = readSegSites(inputFileName, nrHaplotypes, propagationCore.msmc.tTotIntervals);
  foreach(ref s; segsites) {
    if(s.obs.length > 1 || s.obs[0] > 1)
      s.obs = [2];
  }
    
  stderr.writeln("generating Hidden Markov Model");
  return new MSMC_hmm(propagationCore, segsites);

}

MSMC_hmm makeStandardHmm() {
  auto lambdaVec = new double[nrTimeSegments];
  lambdaVec[] = 1.0;
  auto subpopLabels = new size_t[nrHaplotypes];
  auto model = new MSMCmodel(mutationRate, recombinationRate, subpopLabels, lambdaVec, nrTimeSegments, nrTtotSegments);

  stderr.writeln("generating propagation core");
  auto propagationCore = new PropagationCoreFast(model, 1000);

  stderr.writeln("generating Hidden Markov Model");
  return new MSMC_hmm(propagationCore, inputFileName);
}

void decodeWithHmm(MSMC_hmm hmm) {
  hmm.runForward();
  auto forwardState = hmm.propagationCore.newForwardState();
  auto backwardState = hmm.propagationCore.newBackwardState();
  
  double[][] posteriors;
  for(size_t pos = hmm.segsites[$ - 1].pos; pos > hmm.segsites[0].pos && pos <= hmm.segsites[$ - 1].pos; pos -= stride)
  {
    hmm.getForwardState(forwardState, pos);
    hmm.getBackwardState(backwardState, pos);
    auto posterior = forwardState.vecMarginal.dup;
    posterior[] *= backwardState.vecMarginal[];
    auto norm = posterior.reduce!"a+b"();
    posterior[] /= norm;
    posteriors ~= posterior;
  }
  
  foreach_reverse(posterior; posteriors) {
    foreach(p; posterior)
      write(p, " ");
    writeln("");
  }
  
  
}