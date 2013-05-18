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