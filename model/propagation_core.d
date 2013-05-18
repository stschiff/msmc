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
 
module model.propagation_core;
import std.conv;
import model.data;
import model.stateVec;
import model.stateVecAllocator;
import model.msmc_model;

interface PropagationCore {

  void propagateSingleForward(in State_t from, State_t to,
      in SegSite_t from_segsite, in SegSite_t to_segsite) const
  in {
    assert(to_segsite.pos == from_segsite.pos + 1);
  }
  
  void propagateSingleBackward(in State_t to, State_t from,
      in SegSite_t to_segsite, in SegSite_t from_segsite) const
  in {
    assert(to_segsite.pos == from_segsite.pos + 1);
  }
  
  void propagateMultiForward(in State_t from, State_t to,
      in SegSite_t from_segsite, in SegSite_t to_segsite) const
  in {
    assert(to_segsite.pos > from_segsite.pos);
    assert(to_segsite.obs[0] < 2);
  }
  
  void propagateMultiBackward(in State_t to, State_t from,
      in SegSite_t to_segsite, in SegSite_t from_segsite) const
  in {
    assert(to_segsite.pos > from_segsite.pos);
    assert(to_segsite.obs[0] < 2);
  }
  
  @property size_t forwardStateSize() const;
  @property size_t backwardStateSize() const;
  const(MSMCmodel) getMSMC() const;
  State_t newForwardState() const;
  State_t newBackwardState() const;
  State_t newForwardState(StateVecAllocator stateAllocator) const;
  State_t newBackwardState(StateVecAllocator stateAllocator) const;
  void initialState(State_t s) const;
  void getTransitionExpectation(State_t f, State_t b, in SegSite_t to_segsite, double[][] ret) const;
  void setState(State_t s, double x, in SegSite_t segsite) const;
  
  @property size_t maxDistance() const;
  
}