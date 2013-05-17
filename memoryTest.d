import std.stdio;
import propagation_core_fastImpl;
import msmc_model;
import data;
import msmc_hmm;
import stateVecAllocator;

void main() {
  
  auto subpopLabels = [0U, 0, 0, 0];
  auto modelParams = MSMCmodel.withTrivialLambda(0.01, 0.01, subpopLabels, 20, 10);
  
  stderr.writeln("creating propagation core");
  auto propagationCore = new PropagationCoreFast(modelParams, 1000);
  stderr.writeln("done. Press enter");
  stdin.readln();

  stderr.writeln("making hmm");
  auto l = 665621;
  // auto l = 1_000_000;
  // auto stateAllocator = new StateAllocator(l * (propagationCore.forwardStateSize + propagationCore.backwardStateSize));
  SegSite_t[] segsites;
  // size_t[][] vec = new size_t[][l];
  foreach(i; 0 .. l) {
    segsites ~= new SegSite_t(i + 1, [1], 0);
    // vec ~= new size_t[5];
  }
  
  auto hmm = new MSMC_hmm(propagationCore, segsites);
  
  stderr.writeln("done. press enter");
  stdin.readln();
}