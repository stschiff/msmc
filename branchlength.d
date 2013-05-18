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
import std.concurrency;
import std.algorithm;
import std.array;
import std.json;
import std.file;
import std.typecons;
import std.exception;
import std.c.stdlib;
import core.memory;
import model.msmc_hmm;
import model.msmc_model;
import model.triple_index_marginal;
import model.emission_rate;
import model.transition_rate;
import model.time_intervals;
import model.triple_index_marginal;
import model.coalescence_rate;
import model.rate_integrator;
import model.propagation_core_fastImpl;
import model.data;
import utils;
import expectation_step;

double recombinationRate;
double mutationRate;
size_t nrHaplotypes;
string inputFileName;
size_t nrTimeSegments = 40;

void branchlengthMain(string[] args) {
  parseCommandLine(args);
  run();
}

void parseCommandLine(string[] args) {
  if(args.length == 1)
    displayHelpMessageAndExit();
  try {
    readArguments(args);
  }
  catch(Exception e) {
    stderr.writeln("error in parsing command line: ", e.msg);
    displayHelpMessageAndExit();
  }
}
  
void readArguments(string[] args) {
  getopt(args,
      std.getopt.config.caseSensitive,
      "mutationRate|m", &mutationRate,
      "recombinationRate|r", &recombinationRate,
      "nrTimeSegments|n", &nrTimeSegments
  );
  enforce(!isNaN(mutationRate), "need to set mutation rate");
  inputFileName = args[1];
  nrHaplotypes = getNrHaplotypesFromFile(inputFileName);
  if(isNaN(recombinationRate))
    recombinationRate = mutationRate / 4.0;
  stderr.writeln("found ", nrHaplotypes, " Haplotypes in file");
}

static void displayHelpMessageAndExit() {
  stderr.writeln("Usage: msmc branchlength [options] <datafile>
-n, --nrTimeSegments=<int> : nr of time intervals [=40]
-m, --mutationRate=<double> : scaled mutation rate to use
-r, --recombinationRate=<double> : scaled recombination rate to use [default mutationRate / 4]");
  exit(0);
}
  
void run() {
  auto propagationCore = buildPropagationCore();
  auto msmc_hmm = buildHMM(propagationCore);
    
  stderr.writeln("running forward");
  msmc_hmm.runForward();
  
  auto f = File(inputFileName, "r");
  auto forwardState = propagationCore.newForwardState();
  auto backwardState = propagationCore.newBackwardState();

  string[] lines;
  foreach(line; f.byLine()) {
    lines ~= strip(line).idup;
  }
  string[] newLines;
  foreach_reverse(lineIndex; 0 .. lines.length) {
    auto fields = split(lines[lineIndex].dup);
    auto pos = to!int(fields[1]);
    
    msmc_hmm.getForwardState(forwardState, pos);
    msmc_hmm.getBackwardState(backwardState, pos);
    double ttot = 2.0 * propagationCore.msmc.timeIntervals.meanTimeWithLambda(0, 1.0);
    auto max = forwardState.vec[0] * backwardState.vec[0];
    foreach(i; 0 .. nrTimeSegments) {
      auto p = forwardState.vec[i] * backwardState.vec[i];
      if(p > max) {
        max = p;
        // we need the total branch length, so twice the tMRCA with two haplotypes
        ttot = 2.0 * propagationCore.msmc.timeIntervals.meanTimeWithLambda(i, 1.0);
      }
    }
    newLines ~= lines[lineIndex] ~ text(" ", ttot);
  }
  foreach_reverse(line; newLines)
    writeln(line);
}
  
private PropagationCoreFast buildPropagationCore() {
  auto lambdaVec = new double[nrTimeSegments];
  lambdaVec[] = 1.0;
  // the factor 2 is just part of the formula for the mean total branch length.
  auto expectedTtot = 2.0 * TimeIntervals.computeWattersonFactor(nrHaplotypes);
  // the next factor 2 fakes a two haplotype system with the same total branch length (every branch gets half)
  auto boundaries = TimeIntervals.getQuantileBoundaries(nrTimeSegments, expectedTtot / 2.0);
  auto model = new MSMCmodel(mutationRate, recombinationRate, [0UL, 0], lambdaVec, boundaries[0 .. $ - 1], 1);

  stderr.writeln("generating propagation core");
  auto propagationCore = new PropagationCoreFast(model, 1000);
  return propagationCore;
}
  
private MSMC_hmm buildHMM(PropagationCoreFast propagationCore) {
  stderr.writeln("loading file ", inputFileName);
  auto segsites = readSegSites(inputFileName, nrHaplotypes, propagationCore.msmc.tTotIntervals);
  foreach(ref s; segsites) {
    if(s.obs.length > 1 || s.obs[0] > 1)
      s.obs = [2];
  }
    
  stderr.writeln("generating Hidden Markov Model");
  return new MSMC_hmm(propagationCore, segsites);
}
