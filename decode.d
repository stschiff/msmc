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
import std.parallelism;
import std.range;
import model.msmc_hmm;
import model.data;
import model.time_intervals;
import model.msmc_model;
import model.propagation_core_fastImpl;
import branchlength;

double mutationRate, recombinationRate;
size_t nrTimeSegments=40;
size_t nrTtotSegments;
size_t stride=1000;
bool tTot = false;
string inputFileName;
size_t nrHaplotypes;
uint nrThreads;

void main(string[] args) {
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
         "nrThreads", &nrThreads,
         "nrTtotSegments|T", &nrTtotSegments,
         "tTot", &tTot,
         "stride|s", &stride);
  if(nrThreads)
    std.parallelism.defaultPoolThreads(nrThreads);
  enforce(args.length == 2, "need exactly one input file");
  inputFileName = args[1];
  nrHaplotypes = getNrHaplotypesFromFile(inputFileName);
  if(!nrTtotSegments)
      nrTtotSegments = nrTimeSegments;
  enforce(mutationRate > 0, "need positive mutationRate");
  enforce(recombinationRate > 0, "need positive recombinationRate");
}

void displayHelpMessage() {
  stderr.writeln("Usage: msmc decode [options] <inputFile>
Options:
-m, --mutationRate <double>
-r, --recombinationRate <double>
-t, --nrTimeSegments <int> [default=40]
--nrThreads <int> : nr of threads, defaults to nr of CPUs
--tTot : output the total branch length rather than the first coalescent time
-s, --stride <int>: stride width in basepairs [default=1000]");
}

void run() {
  auto hmm = tTot ? makeTtotHmm() : makeStandardHmm();
  decodeWithHmm(hmm);
}

MSMC_hmm makeTtotHmm() {
  auto lambdaVec = new double[nrTtotSegments];
  lambdaVec[] = 1.0;
  auto expectedTtot = 2.0;
  auto boundaries = TimeIntervals.getQuantileBoundaries(nrTtotSegments, expectedTtot / 2.0);
  auto model = new MSMCmodel(mutationRate, recombinationRate, [0UL, 0], lambdaVec, boundaries[0 .. $ - 1], 1, false);

  stderr.writeln("generating propagation core");
  auto propagationCore = new PropagationCoreFast(model, 1000);

  stderr.writeln("loading file ", inputFileName);
  auto segsites = readSegSites(inputFileName, false, [], false);
  auto allele_order = canonicalAlleleOrder(nrHaplotypes);
  foreach(ref s; segsites) {
    if(s.obs[0] > 1) {
      auto count_1 = count(allele_order[s.obs[0] - 1], '1');
      if(count_1 == 1)
        s.obs = [2];
      else
        s.obs = [1];
    }
  }
  stderr.writeln("generating Hidden Markov Model");
  return new MSMC_hmm(propagationCore, segsites);

}

MSMC_hmm makeStandardHmm() {
  auto lambdaVec = new double[nrTimeSegments];
  lambdaVec[] = 1.0;
  auto subpopLabels = new size_t[nrHaplotypes];
  auto model = new MSMCmodel(mutationRate, recombinationRate, subpopLabels, lambdaVec, nrTimeSegments, nrTtotSegments, false);

  stderr.writeln("generating propagation core");
  auto propagationCore = new PropagationCoreFast(model, 1000);
  
  auto segsites = readSegSites(inputFileName, false, [], false);
  if(nrHaplotypes > 2)
    estimateTotalBranchlengths(segsites, model);

  stderr.writeln("generating Hidden Markov Model");
  return new MSMC_hmm(propagationCore, segsites);
}

void decodeWithHmm(MSMC_hmm hmm) {
  hmm.runForward();
  auto forwardState = hmm.propagationCore.newForwardState();
  auto backwardState = hmm.propagationCore.newBackwardState();
  
  double[][] posteriors;
  double [] positions;
  for(size_t pos = hmm.segsites[$ - 1].pos; pos > hmm.segsites[0].pos && pos <= hmm.segsites[$ - 1].pos; pos -= stride)
  {
    hmm.getForwardState(forwardState, pos);
    hmm.getBackwardState(backwardState, pos);
    auto posterior = forwardState.vecMarginal.dup;
    posterior[] *= backwardState.vecMarginal[];
    auto norm = posterior.reduce!"a+b"();
    posterior[] /= norm;
    posteriors ~= posterior;
    positions ~= pos;
  }
  
  auto timeBoundaries = hmm.propagationCore.getMSMC().timeIntervals.boundaries.dup;
  timeBoundaries[0] = 0.0;
  write("#time_boundaries:\t");
  foreach(t; timeBoundaries[0 .. $ - 1])
      write(t, "\t");
  writeln("");
  
  foreach_reverse(e; zip(positions, posteriors)) {
    write(e[0], "\t");
    foreach(p; e[1])
      write(p, "\t");
    writeln("");
  }
  
  
}