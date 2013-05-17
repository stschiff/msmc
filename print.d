import std.stdio;
import std.json;
import std.file;
import std.exception;
import std.getopt;
import std.c.stdlib;
import utils;

JSONValue json;
enum what_t {popSize, mutationRate, recombinationRate, lambda, crossLambda, likelihood};
what_t what;
bool unscaled;
bool years;
double mu = 1.25e-8;
char delimiter = '\t';

void printMain(string[] args) {
  try {
    readArgs(args);
  }
  catch (Exception e) {
    stderr.writeln(e.msg);
    printHelpMessage();
    exit(0);
  }
  final switch(what) {
    case what_t.popSize:
    runPopSize();
    break;
    case what_t.mutationRate:
    runMutationRate();
    break;
    case what_t.recombinationRate:
    runRecombinationRate();
    break;
    case what_t.lambda:
    runLambda();
    break;
    case what_t.crossLambda:
    runCrossLambda();
    break;
    case what_t.likelihood:
    runLikelihood();
    break;
  }
}

void readArgs(string[] args) {
  getopt(args, "what|w", &what, "unscaled|s", &unscaled, "years|y", &years, "delimiter|d", &delimiter);
  enforce(args.length == 2, "need exactly one filename");
  json = parseJSON(readText(args[1]));
}

void printHelpMessage() {
  stderr.writeln("./jsonToTable [Options] filename
    Options:
    --what, -w <what>: one of [popSize, mutationRate, recombinationRate, lambda, crossLambda, likelihood]. Default [popSize]
    --unscaled, -s
    --years, -y print time in years, not generations (only applied if --unscaled is set)
    --delimiter, -d");
}

void runPopSize() {
  auto times = json["results"].array[$ - 1]["updatedParams"]["timeIntervals"].fromJSON!(double[])();
  auto lambdaVec = json["results"].array[$ - 1]["updatedParams"]["lambdaVec"].fromJSON!(double[])();
  enforce(lambdaVec.length == times.length, "cross-population run, use lambda instead");
  times[0] = times[1] / 4.0;
  foreach(i; 0 .. lambdaVec.length) {
    if(unscaled) {
      auto t = 2.0 * getN() * times[i];
      if(years)
        t *= 30;
      writeln(t, delimiter, getN() / lambdaVec[i]);
    }
    else
      writeln(times[i], delimiter, 1.0 / lambdaVec[i]);
  }
}

void runCrossLambda() {
  auto times = json["results"].array[$ - 1]["updatedParams"]["timeIntervals"].fromJSON!(double[])();
  auto lambdaVec = json["results"].array[$ - 1]["updatedParams"]["lambdaVec"].fromJSON!(double[])();
  enforce(lambdaVec.length == 3 * times.length, "need cross-population run with 2 populations");
  times[0] = times[1] / 4.0;
  foreach(i; 0 .. times.length) {
    if(unscaled) {
      auto t = 2.0 * getN() * times[i];
      if(years)
        t *= 30;
      writeln(t, delimiter, 2.0 * lambdaVec[3 * i + 1] / (lambdaVec[3 * i] + lambdaVec[3 * i + 2]));
    }
    else
      writeln(times[i], delimiter, 2.0 * lambdaVec[3 * i + 1] / (lambdaVec[3 * i] + lambdaVec[3 * i + 2]));
  }
}

void runLambda() {
  auto times = json["results"].array[$ - 1]["updatedParams"]["timeIntervals"].fromJSON!(double[])();
  auto lambdaVec = json["results"].array[$ - 1]["updatedParams"]["lambdaVec"].fromJSON!(double[])();
  times[0] = times[1] / 4.0;
  auto nrSubpopPairs = lambdaVec.length / times.length;
  foreach(i; 0 .. times.length) {
    if(unscaled) {
      auto t = 2.0 * getN() * times[i];
      if(years)
        t *= 30;
      write(t);
    }
    else
      write(times[i]);
    foreach(subpopPair; 0 .. nrSubpopPairs) {
      write(delimiter, lambdaVec[nrSubpopPairs * i + subpopPair]);
    }
    writeln("");
  }
}

void runMutationRate() {
  auto scaledMut = json["results"].array[$ - 1]["updatedParams"]["mutationRate"].fromJSON!double();
  writeln(scaledMut);
}

void runRecombinationRate() {
  foreach(res; json["results"].array) {
    auto rho = res["updatedParams"]["recombinationRate"].fromJSON!double();
    if(unscaled)
      rho /= (2.0 * getN());
    writeln(rho);
  }
}

void runLikelihood() {
  foreach(res; json["results"].array) {
    auto likelihood = res["logLikelihood"].fromJSON!double();
    writeln(likelihood);
  }
}

double getN() {
  auto scaledMu = json["results"].array[$ - 1]["updatedParams"]["mutationRate"].fromJSON!double();
  return scaledMu / (2.0 * mu);
}