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
  
