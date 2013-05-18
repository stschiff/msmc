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

module model.stateVec;
import std.conv;
import std.algorithm;
import model.gsl_matrix_vector;
import model.stateVecAllocator;

class State_t {
  double[] vec, vecMarginal, vecMarginalEmission;
  
  this(size_t nrS, size_t nrM, size_t nrME) {
    vec = new double[nrS];
    vecMarginal = new double[nrM];
    vecMarginalEmission = new double[nrME];
  }
  
  this(size_t nrS, size_t nrM, size_t nrME, StateVecAllocator stateAllocator) {
    vec = stateAllocator.allocate(nrS);
    vecMarginal = stateAllocator.allocate(nrM);
    vecMarginalEmission = stateAllocator.allocate(nrME);
  }
  
  @property double norm() const {
    if(vecMarginal.length > 0)
      return reduce!"a+b"(vecMarginal);
    else
       return reduce!"a+b"(vec);
  }
  
  void scale(double x) {
    vec[] *= x;
    vecMarginal[] *= x;
    vecMarginalEmission[] *= x;
  }
  
  void copy_into(State_t dest) {
    dest.vec[] = vec[];
    dest.vecMarginal[] = vecMarginal[];
    dest.vecMarginalEmission[] = vecMarginalEmission[];
  }
  
  override string toString() const {
    return text(vec);
  }

  void setZero() {
    vec[] = 0.0;
    vecMarginal[] = 0.0;
    vecMarginalEmission[] = 0.0;
  }

}
