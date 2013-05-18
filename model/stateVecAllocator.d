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
 
class StateVecAllocator {
  double[] vec;
  size_t bytesUsed;
  
  this(size_t reservedSize) {
    vec = new double[reservedSize];
    bytesUsed = 0;
  }
  
  double[] allocate(size_t size) {
    auto start = bytesUsed;
    auto end = start + size;
    if(vec.length < end)
      throw new Exception("State Allocator full");
    bytesUsed += size;
    return vec[start .. end];
  }
}

unittest {
  import std.stdio;
  import std.exception;
  stderr.writeln("testing StateAllocator");
  
  auto stateAllocator =  new StateVecAllocator(100);
  
  auto vec = stateAllocator.allocate(90);
  assertThrown(stateAllocator.allocate(20));
}