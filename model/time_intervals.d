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
 
import std.exception;
import std.stdio;
import std.conv;
import std.math;
import std.algorithm;
import triple_index_marginal;

class TimeIntervals {
  const double[] boundaries;
  alias BoundaryFunction_t = double function(size_t, size_t, double);
  
  static double computeLiAndDurbinBoundary(size_t i, size_t nrTimeSegments, double factor) {
    auto tMax = 15.0;
    auto boundary = factor * (0.1 * (exp(i / cast(double)nrTimeSegments * log(1.0 + 10.0 * tMax)) - 1.0));
    if(i == nrTimeSegments)
      boundary = double.infinity;
    return boundary;
  }
  
  static double computeQuantileBoundary(size_t i, size_t nrTimeSegments, double factor) {
    return -factor * log(1.0 - cast(double)i / cast(double)nrTimeSegments);
  }
  
  static double[] getBoundaries(BoundaryFunction_t boundaryFunction, size_t nrTimeSegments, double factor) {
    enforce(factor > 0.0);
    enforce(nrTimeSegments > 0);
    auto boundaries = new double[nrTimeSegments + 1];
    foreach(i; 0 .. nrTimeSegments + 1)
      boundaries[i] = boundaryFunction(i, nrTimeSegments, factor);
    return boundaries;
  }
  
  static double[] getQuantileBoundaries(size_t nrTimeSegments, double factor) {
    return getBoundaries(&computeQuantileBoundary, nrTimeSegments, factor);
  }

  static double[] getLiAndDurbinBoundaries(size_t nrTimeSegments, double factor) {
    return getBoundaries(&computeLiAndDurbinBoundary, nrTimeSegments, factor);
  }
  
  static TimeIntervals standardIntervals(size_t nrTimeSegments, size_t nrHaplotypes) {
    auto mOver2 = nrHaplotypes * (nrHaplotypes - 1) / 2;
    auto boundaries = getBoundaries(&computeQuantileBoundary, nrTimeSegments, 1.0 / mOver2);
    return new TimeIntervals(boundaries);
  }

  static TimeIntervals standardTotalBranchlengthIntervals(size_t nrTimeSegments, size_t nrHaplotypes) {
    auto expectedTtot = 2.0 * computeWattersonFactor(nrHaplotypes);
    auto boundaries = getBoundaries(&computeQuantileBoundary, nrTimeSegments, expectedTtot);
    return new TimeIntervals(boundaries);
  }
  
  static double computeWattersonFactor(size_t nrHaplotypes) {
    auto wattersonFactor = 0.0;
    foreach(i; 1 .. nrHaplotypes)
      wattersonFactor += 1.0 / i;
    return wattersonFactor;
  }
  
  this(in double[] boundaries)
  {
    enforce(validBoundaries(boundaries), text("invalid boundaries: ", boundaries));
    this.boundaries = boundaries.dup;
  };
  
  private static bool validBoundaries(in double[] tb) {
    if(tb[0] != 0.0 || tb[$ - 1] != double.infinity)
      return false;
    foreach(i; 0 .. tb.length - 1)
      if(tb[i] >= tb[i + 1])
        return false;
    return true;
  }

  double leftBoundary(size_t i) const {
    return boundaries[i];
  }
  
  double rightBoundary(size_t i) const {
    return boundaries[i + 1];
  }
  
  double delta(size_t i) const {
    return boundaries[i + 1] - boundaries[i];
  }
  
  size_t findIntervalForTime(double t) const
    in {
      assert(t >= 0.0);
    }
    out(result) {
      assert(result < nrIntervals());
    }
  body {
    auto index = countUntil!"a>b"(boundaries[1 .. $], t);
    if(index < 0)
      index = nrIntervals - 1;
    return index;
  }
  
  size_t roundToFullInterval(double t) const
    in {
      assert(t >= 0.0);
    }
  body {
    auto index = findIntervalForTime(t);
    if(t  - boundaries[index] < boundaries[index + 1] - t)
      return index;
    else
      return index + 1;
  }
  
  @property size_t nrIntervals() const {
    return cast(size_t)boundaries.length - 1;
  }
  
  double meanTimeWithLambda(size_t i, double lambda) const
    out(result) {
      assert(result > leftBoundary(i) && result < rightBoundary(i));
    }
  body {
    double time;
    if(lambda < 0.001) {
      return meanTimeForLowRate(i, lambda);
    }
    else {
      return meanTimeForNormalRate(i, lambda);
    }
  }
  
  private double meanTimeForLowRate(size_t i, double lambda) const {
    if(i < nrIntervals() - 1)
      return (leftBoundary(i) + rightBoundary(i)) / 2.0;
    else
      return leftBoundary(i) + 1.0 / lambda;
  }

  private double meanTimeForNormalRate(size_t i, double lambda) const {
    auto time = 1.0 + lambda * leftBoundary(i);
    if(i < nrIntervals() - 1)
      time -= exp(-lambda * delta(i)) * (1.0 + lambda * rightBoundary(i));
    time /= (1.0 - exp(-delta(i) * lambda)) * lambda;
    return time;
  }
  
  double meanTime(size_t i) const {
    return meanTimeForLowRate(i, 1.0);
  }
}

double[] convertPiecewiseFunctions(in TimeIntervals timeIntervals, in double[2][] functionPoints) {
  auto newIntervalValues = new double[timeIntervals.nrIntervals];
  
  foreach(i; 0 .. functionPoints.length) {
    auto index = timeIntervals.roundToFullInterval(functionPoints[i][0]);
    newIntervalValues[index .. $] = functionPoints[i][1];
  }
  return newIntervalValues;
}

unittest {
  writeln("test convertPiecewiseFunctions");
  auto boundaries = [0.0, 1.0, 4.0, 10.0, 20.0, 30.0, double.infinity];
  auto timeIntervals = new TimeIntervals(boundaries);
  double[2][] functionPoints = [[0.0, 3.0], [0.1, 6.0], [12.0, 2.0], [200.0, 6.7]];
  auto newIntervalValues = convertPiecewiseFunctions(timeIntervals, functionPoints);
  assert(newIntervalValues[0] == 6.0);
  assert(newIntervalValues[1] == 6.0);
  assert(newIntervalValues[2] == 6.0);
  assert(newIntervalValues[3] == 2.0);
  assert(newIntervalValues[4] == 2.0);
  assert(newIntervalValues[5] == 6.7);
}

unittest {
  writeln("test timeInterval boundaries and delta");
  auto t = TimeIntervals.standardIntervals(10, 4);
  
  assert(t.leftBoundary(0) == 0);
  assert(t.rightBoundary(9) == double.infinity);
  assert(t.delta(4) == t.rightBoundary(4) - t.leftBoundary(4));
}

unittest {
  writeln("test timeIntervals.findIntervalForTime");
  auto t = TimeIntervals.standardIntervals(10, 4);
  auto t1 = t.leftBoundary(4);
  auto t2 = t.rightBoundary(4);
  auto middle = (t1 + t2) / 2.0;
  
  assert(t.findIntervalForTime(middle) == 4);
  assert(t.findIntervalForTime(t1) == 4);
  assert(t.findIntervalForTime(t2) == 5);
  assert(t.findIntervalForTime(0.0) == 0);
  assert(t.findIntervalForTime(double.infinity) == 9);
}

unittest {
  writeln("test timeIntervals.roundToFullInterval");
  auto t = TimeIntervals.standardIntervals(10, 4);
  auto t1 = t.leftBoundary(4);
  auto t2 = t.rightBoundary(4);
  
  assert(t.roundToFullInterval(t1-1e-8) == 4);
  assert(t.roundToFullInterval(t1+1e-8) == 4);
  assert(t.roundToFullInterval(t2-1e-8) == 5);
}

unittest {
  writeln("test timeIntervals.validBoundaries");
  assert(TimeIntervals.validBoundaries([0.0, 1.0, 4.0, double.infinity]));
  assert(!TimeIntervals.validBoundaries([1.0, double.infinity]));
  assert(!TimeIntervals.validBoundaries([0.0, 3.0, 2.0, double.infinity]));
}

unittest {
  writeln("test timeIntervals.meanTime");
  import coalescence_rate;
  auto meanTimesTh = [0.02053010050230513,0.10097875411723184,0.2850866082677072,0.8012622730318943];
  auto subpop_labels = [0UL, 0, 1, 1];
  auto lambda_subpop_rates = [1, 0.1, 1, 2, 0.5, 4, 2, 1.0, 4, 1., 2, 1];
  auto T = 4UL;
  auto marginalIndex = new MarginalTripleIndex(T, subpop_labels);
  auto coal = new PiecewiseConstantCoalescenceRate(marginalIndex, lambda_subpop_rates);
  auto boundaries = TimeIntervals.getLiAndDurbinBoundaries(T, 1.0 / 6.0);
  auto timeIntervals = new TimeIntervals(boundaries);
  foreach(i; 0 .. T)
    assert(approxEqual(timeIntervals.meanTimeWithLambda(i, coal.getTotalMarginalLambda(i)), meanTimesTh[i], 1.0e-8, 0.0));
}

unittest {
  writeln("test one time interval");
  auto t = TimeIntervals.standardIntervals(1, 2);
  assert(t.leftBoundary(0) == 0);
  assert(t.rightBoundary(0) == double.infinity);
  assert(t.delta(0) == double.infinity);
  assert(t.findIntervalForTime(0.0) == 0);
  assert(t.findIntervalForTime(1.0) == 0);
  auto tMean = t.meanTimeWithLambda(0, 1.0);
  assert(tMean > 0 && tMean < double.infinity);
}

unittest {
  writeln("test boundary infinity");
  assert(TimeIntervals.computeQuantileBoundary(10, 10, 1.0) == double.infinity);
  assert(TimeIntervals.computeLiAndDurbinBoundary(10, 10, 1.0) == double.infinity);
}