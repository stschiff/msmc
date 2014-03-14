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
import std.string;
import std.array;
import model.msmc_model;
import maximization_step;

void main(string[] args) {
  auto mu = 0.000336401;
  auto rho = 8.41002e-05;
  auto timeSegments = [1UL, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
  auto T = timeSegments.reduce!"a+b"();
  auto subpopLabels = [0UL, 0, 1, 1];
  
  auto msmc = MSMCmodel.withTrivialLambda(mu, rho, subpopLabels, T, 30, false);

  string filename = args[1];
  
  auto f = File(filename, "r");
  auto line = f.readln();
  auto eVec = line.strip.split.map!"to!double(a)"().array();
  double[][] eMat;
  foreach(l; f.byLine)
    eMat ~= l.strip.split.map!"to!double(a)"().array();
  
  // stderr.writeln(eVec);
  // stderr.writeln(eMat);
  
  auto newParams = getMaximization(eVec, eMat, msmc, timeSegments, false, true);
  
  foreach(i; 0 .. T) {
    writefln("1\t%s\t%s", i + 1, newParams.lambdaVec[i]);
  }
  
}
