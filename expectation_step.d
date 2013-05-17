import std.typecons;
import std.concurrency;
import std.stdio;
import std.string;
import std.algorithm;
import std.exception;
import core.memory;
import propagation_core;
import propagation_core_fastImpl;
import propagation_core_naiveImpl;
import msmc_model;
import msmc_hmm;

class ExpectationStep {
  
  alias immutable(immutable(double)[])[] immutableMatrix;
  
  private const MSMCmodel msmc;
  private immutable PropagationCore propagationCore;
  private const string[] fileNames;
  private size_t hmmStrideWidth;
  private double[][] expectationResult;
  private double logLikelihood;
  private size_t nrSubThreads;
  
  this(size_t nrSubThreads, in string[] fileNames, MSMCmodel msmc, size_t hmmStrideWidth,
       size_t maxDistance=1000, bool naiveImplementation=false)
  {
    enforce(nrSubThreads > 0, "need at least one thread");
    enforce(fileNames.length > 0, "need at least one file");
    enforce(hmmStrideWidth > 0, "need positive hmm stride width");
    this.msmc = msmc;
    this.fileNames = fileNames.dup;
    this.hmmStrideWidth = hmmStrideWidth;
    this.nrSubThreads = nrSubThreads;
    stderr.writeln("building propagation core");
    if(naiveImplementation)
      propagationCore = new immutable(PropagationCoreNaive)(msmc, maxDistance);
    else
      propagationCore = new immutable(PropagationCoreFast)(msmc, maxDistance);
    stderr.writeln("propagation core ready");
    initResults();
  }
  
  private void initResults() {
    expectationResult = new double[][](msmc.nrMarginals, msmc.nrMarginals);
    foreach(ref row; expectationResult)
      row[] = 0.0;
    logLikelihood = 0.0;
  }
  
  void run() {
    if(nrSubThreads == 1)
      runSingleThreaded();
    else
      runMultiThreaded();
    cleanup();
  }
  
  private void runSingleThreaded() {
    foreach(filename; fileNames) {
      auto result = singleChromosomeExpectation(filename, hmmStrideWidth, propagationCore);
      foreach(au; 0 .. msmc.nrMarginals)
        expectationResult[au][] += result[0][au][];
      logLikelihood += result[1];
    }
  }
  
  private void runMultiThreaded() {
    auto index = 0;
    foreach(i; 0 .. min(nrSubThreads, fileNames.length)) {
      spawn(&singleChromosomeExpectationThread, thisTid, fileNames[index], propagationCore, hmmStrideWidth);
      index += 1;
    }
    
    for(int resultsReceived = 0; resultsReceived < fileNames.length; ) {
      receive(
        (Tuple!(immutableMatrix, double) msg) {
          stderr.writefln("received expectation result");
          auto transitions = msg[0];
          logLikelihood += msg[1];
          foreach(au; 0 .. msmc.nrMarginals)
            expectationResult[au][] += transitions[au][];
          resultsReceived += 1;
          if(index < fileNames.length) {
            spawn(&singleChromosomeExpectationThread, thisTid, fileNames[index], propagationCore, hmmStrideWidth);
            index += 1;
          }
        },
        (LinkTerminated e) {
          writeln("caught Exception", e.msg);
        }
      );
    }
  }
  
  private static Tuple!(immutableMatrix, double)
  singleChromosomeExpectation(string filename, size_t hmmStrideWidth, in PropagationCore propagationCore)
  {
    auto logTag = format("Subthread for file %s: ", filename);
    auto msmc_hmm =  new MSMC_hmm(propagationCore, filename, logTag);

    stderr.writeln(logTag ~ "running forward");
    msmc_hmm.runForward();
    stderr.writeln(logTag ~ "running backward");
    auto exp = msmc_hmm.runBackward(hmmStrideWidth);
    auto logL = msmc_hmm.logLikelihood();
    stderr.writeln(logTag, "likelihood: ", logL);
    msmc_hmm.destroy();
    return tuple(exp, logL);
  }

  private static void singleChromosomeExpectationThread(Tid caller, string fileName,
      immutable PropagationCore propagationCore, size_t hmmStrideWidth)
  {
    auto result = singleChromosomeExpectation(fileName, hmmStrideWidth, propagationCore);
    caller.send(result[0], result[1]);
  }
  
  private void cleanup() {
    propagationCore.destroy();
    GC.collect();
    GC.minimize();
  }
  
  double[][] getResult() {
    return expectationResult;
  }
  
  double getLogLikelihood() {
    return logLikelihood;
  }

}

unittest {
  writeln("test expectation step");
  auto lambdaVec = new double[12];
  lambdaVec[] = 1.0;
  auto msmc = new MSMCmodel(0.01, 0.001, [0U, 0, 1, 1], lambdaVec, 4, 4);
  auto fileNames = ["model/hmm_testData_tMRCA.txt", "model/hmm_testData_tMRCA.txt"];
  auto hmmStrideWidth = 100U;
  
  auto expectationStep = new ExpectationStep(1U, fileNames, msmc, hmmStrideWidth, 100U);
  expectationStep.run();
  auto resultSingleThreaded = expectationStep.getResult();
  expectationStep = new ExpectationStep(2U, fileNames, msmc, hmmStrideWidth, 100U);
  expectationStep.run();
  auto resultMultiThreaded = expectationStep.getResult();
   
  assert(resultSingleThreaded == resultMultiThreaded);
  foreach(row; resultSingleThreaded)
    foreach(val; row)
      assert(val >= 0.0);
}
