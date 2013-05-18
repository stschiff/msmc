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
import std.algorithm;
import std.array;
import std.conv;
import std.c.stdlib;
import inference;
import expectation;
import maximization;
import stats;
import branchlength;
import decode;
import print;

void main(string[] args) {
  if(args.length < 2) {
    exitWithHelpMessage();
  }
  switch(args[1]) {
    case "inference":
    inferenceMain(args[1..$]);
    break;
    case "branchlength":
    branchlengthMain(args[1..$]);
    break;
    case "stats":
    StatsApplication.runWithArgs(args[1..$]);
    break;
    case "expectation":
    ExpectationApplication.runWithArgs(args[1..$]);
    break;
    case "maximization":
    MaximizationApplication.runWithArgs(args[1..$]);
    break;
    case "decode":
    decodeMain(args[1..$]);
    break;
    case "print":
    printMain(args[1..$]);
    break;
    default:
    writeln("subprogram ", args[1], " not recognized");
    exitWithHelpMessage();
    break;
  }
}

void exitWithHelpMessage() {
  writeln("Usage:
msmc inference - infer population size history and migration patterns from multiple haplotypes
msmc expectation - carry out only one expectation step (forward-backward)
msmc maximization - carry out only one maximization step
msmc branchlengh - annotated a datafile with estimates of the local total branch length
msmc stats - measure diversity stats in multiple datafiles
msmc decode - output the posterior probability of each state locally
msmc print - print various results from inference json output");
  exit(0);
}
  
