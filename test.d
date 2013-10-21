import std.stdio;
import model.msmc_model;
import model.propagation_core_fastImpl;

void main() {
  auto mu = 0.0003493;
  auto rho = 8.7325e-05;
  auto subpopLabels = [0UL, 0, 0, 0];
  auto nrS = 120;
  auto nrObs = 8;

  auto lambda0 = [1.0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
  auto lambda1 = [2.72772,4.82693,4.6929,4.19957,3.58957,3.075,2.65873,2.29462,1.98028,1.71458,1.4953,1.3175,1.164,1.02722,0.905408,0.793915,0.684908,0.573413,0.455717,0.000143077];
  
  auto params0 = new MSMCmodel(mu, rho, subpopLabels, lambda0, 20, 10);
  auto params1 = new MSMCmodel(mu, rho, subpopLabels, lambda1, 20, 10);
  
  auto prop0 = new PropagationCoreFast(params0, 1000);
  auto prop1 = new PropagationCoreFast(params1, 1000);
  
  foreach(i_ttot; 0 .. 10) {
    foreach(aij; 0 .. nrS) {
      auto sum0 = 0.0;
      auto sum1 = 0.0;
      foreach(obs; 1 .. nrObs + 1) {
        sum0 += prop0.emissionProbs[i_ttot][obs][aij];
        sum1 += prop1.emissionProbs[i_ttot][obs][aij];
      }
      writefln("%s\t%s\t%s\t%s", i_ttot, aij, sum0, sum1);
    }
  }

  // foreach(bv; 0 .. 20) {
  //   auto sum0 = 0.0;
  //   auto sum1 = 0.0;
  //   foreach(au; 0 .. 20) {
  //     sum0 += 6 * prop0.transitionMatrixQ2[au][bv];
  //     sum1 += 6 * prop1.transitionMatrixQ2[au][bv];
  //   }
  //   sum0 += prop0.transitionMatrixQ1[bv];
  //   sum1 += prop1.transitionMatrixQ1[bv];
  //   writefln("%s\t%s\t%s", bv, sum0, sum1);
  // }
  
}