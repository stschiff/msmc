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
 
#!/usr/bin/env rdmd

import std.stdio;
import std.algorithm;
import std.string;
import std.array;
import model.msmc_model;
import maximization_step;

void main() {
  auto mu = 0.0003578;
  auto rho = 0.0001;
  auto timeSegments = [1UL, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2];
  auto T = timeSegments.reduce!"a+b"();
  auto subpopLabels = [0UL, 0, 0, 0, 0, 0];
  
  auto msmc = MSMCmodel.withTrivialLambda(mu, rho, subpopLabels, T, 20);

  string filename = "/Users/ss27/Data/MSMC_simulations_new/M6_bottleneck_fine/simState_transitions.txt";
  
  auto f = File(filename, "r");
  
  auto eVec = new double[T];
  auto eMat = new double[][](T, T);
  eVec[] = 0.0;
  foreach(ref r; eMat)
    r[] = 0.0;
  
  double[][] mat;
  foreach(line; f.byLine) {
    mat ~= line.strip.split.map!"to!double(a)"().array();
  }
  
  foreach(aij; 0 .. msmc.nrStates) {
    auto au = msmc.marginalIndex.getMarginalIndexFromIndex(aij);

    foreach(bkl; 0 .. msmc.nrStates) {
      auto bv = msmc.marginalIndex.getMarginalIndexFromIndex(bkl);
      if(aij != bkl) {
        eMat[au][bv] += mat[aij][bkl];
      }
      else
        eVec[au] += mat[aij][bkl];
    }
  }
  auto newParams = getMaximization(eVec, eMat, msmc, timeSegments, false, true);
  
  foreach(i; 0 .. T) {
    writefln("1\t%s\t%s", i + 1, newParams.lambdaVec[i]);
  }
  
}
