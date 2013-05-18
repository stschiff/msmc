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
import model.propagation_core_fastImpl;
import model.msmc_model;
import model.data;
import model.msmc_hmm;
import model.stateVecAllocator;

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