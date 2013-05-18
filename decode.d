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
import std.getopt;
import std.exception;
import std.c.stdlib;
import std.algorithm;
import model.msmc_hmm;
import model.data;
import model.time_intervals;
import model.msmc_model;
import model.propagation_core_fastImpl;

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