import time_intervals;
import coalescence_rate;
import triple_index_marginal;
import std.math;
import std.stdio;

class CoalescenceRateIntegrator {
  const TimeIntervals timeIntervals;
  const PiecewiseConstantCoalescenceRate coal;
  
  this(in TimeIntervals timeIntervals, in PiecewiseConstantCoalescenceRate coal) {
    this.timeIntervals = timeIntervals;
    this.coal = coal;
  }
    
  private double LIntegral(double delegate(size_t) lambdaFunc, double from, double to) const
  {
    auto fromIndex = timeIntervals.findIntervalForTime(from);
    auto toIndex = timeIntervals.findIntervalForTime(to);
    double sum = 0.0;
    if(fromIndex == toIndex) {
      return exp(-(to - from) * lambdaFunc(toIndex));
    }
    foreach(kappa; fromIndex + 1 .. toIndex) {
      sum += lambdaFunc(kappa) * timeIntervals.delta(kappa);
    }
    
    double ret = exp(-(timeIntervals.rightBoundary(fromIndex) - from) *
               lambdaFunc(fromIndex) - sum - (to - timeIntervals.leftBoundary(toIndex)) * 
               lambdaFunc(toIndex));
    
    return ret;
  }
  
  double integrateSemiMarginalLambda(size_t branch, double from, double to) const
  {
    double lambdaFunc(size_t a) {
      return coal.getSemiMarginalLambda(a, branch);
    }
    return LIntegral(&lambdaFunc, from, to);
  }

  double integrateTotalMarginalLambda(double from, double to) const {
    return LIntegral(&coal.getTotalMarginalLambda, from, to);
  }
  
}