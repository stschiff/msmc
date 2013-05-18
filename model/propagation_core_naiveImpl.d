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
 
// module msmc.propagation_core_naiveImpl;
import std.stdio;
import std.algorithm;
import std.conv;
import std.string;
import std.exception;
import data;
import gsl_matrix_vector;
import propagation_core;
import msmc_model;
import stateVec;
import stateVecAllocator;

class PropagationCoreNaive : PropagationCore {
  
  gsl_matrix*[][] forwardPropagators, backwardPropagators;
  gsl_matrix*[] forwardPropagatorsMissing, backwardPropagatorsMissing;
  
  gsl_vector*[][] emissionProbs; // first index: i_ttot, second index: obs
  gsl_matrix* transitionMatrix;
  
  const MSMCmodel msmc;
  
  this(in MSMCmodel msmc, size_t maxDistance) {
    enforce(maxDistance > 0);
    this.msmc = msmc;
      
    auto allele_order = canonicalAlleleOrder(msmc.nrHaplotypes);
    forwardPropagators = new gsl_matrix*[][](msmc.nrTtotIntervals, maxDistance);
    backwardPropagators = new gsl_matrix*[][](msmc.nrTtotIntervals, maxDistance);
    emissionProbs = new gsl_vector*[][](msmc.nrTtotIntervals, allele_order.length + 1);
      
    foreach(tt; 0 .. msmc.nrTtotIntervals) {
      foreach(i; 0 .. allele_order.length + 1) {
        emissionProbs[tt][i] = gsl_vector_alloc(msmc.nrStates);
        foreach(aij; 0 .. msmc.nrStates) {
          if(i == 0)
            gsl_vector_set(emissionProbs[tt][i], aij, 1.0); // missing data
          else
            gsl_vector_set(emissionProbs[tt][i], aij, msmc.emissionProb(allele_order[i - 1], aij, tt));
        }
      }
    }
    
    transitionMatrix = gsl_matrix_alloc(msmc.nrStates, msmc.nrStates);
    foreach(aij; 0 .. msmc.nrStates) {
      foreach(bkl; 0 .. msmc.nrStates) {
        gsl_matrix_set(transitionMatrix, aij, bkl, msmc.transitionRate.transitionProbability(aij, bkl));
      }
    }
      
    foreach(i; 0 .. msmc.nrTtotIntervals) {
      // stderr.writeln("computing propagators, i_ttot=", i);
      foreach(dist; 0 .. maxDistance) {
        forwardPropagators[i][dist] = gsl_matrix_alloc(msmc.nrStates, msmc.nrStates);
        backwardPropagators[i][dist] = gsl_matrix_alloc(msmc.nrStates, msmc.nrStates);
      }
      computeForwardPropagators(forwardPropagators[i], false, i, maxDistance);
      computeBackwardPropagators(backwardPropagators[i], false, i, maxDistance);
    }
    // stderr.writeln("computing propagators for missing data");
    forwardPropagatorsMissing = new gsl_matrix*[maxDistance];
    backwardPropagatorsMissing = new gsl_matrix*[maxDistance];
    foreach(dist; 0 .. maxDistance) {
      forwardPropagatorsMissing[dist] = gsl_matrix_alloc(msmc.nrStates, msmc.nrStates);
      backwardPropagatorsMissing[dist] = gsl_matrix_alloc(msmc.nrStates, msmc.nrStates);
    }
    computeForwardPropagators(forwardPropagatorsMissing, true, 0, maxDistance);
    computeBackwardPropagators(backwardPropagatorsMissing, true, 0, maxDistance);
    
  }
  
  ~this() {
    foreach(i; 0 .. forwardPropagators.length) {
      foreach(dist; 0 .. forwardPropagators[i].length) {
        gsl_matrix_free(forwardPropagators[i][dist]);
        gsl_matrix_free(backwardPropagators[i][dist]);
      }
    }
    foreach(dist; 0 .. forwardPropagatorsMissing.length) {
      gsl_matrix_free(forwardPropagatorsMissing[dist]);
      gsl_matrix_free(backwardPropagatorsMissing[dist]);
    }
    foreach(tt; 0 .. emissionProbs.length)
      foreach(i; 0 .. emissionProbs[tt].length)
        gsl_vector_free(emissionProbs[tt][i]);
    gsl_matrix_free(transitionMatrix);
  }
  
  private void computeForwardPropagators(gsl_matrix*[] ret,
      bool missing_data, size_t i_Ttot, size_t maxDistance) const
  {
    foreach(aij; 0 .. msmc.nrStates) {
      double e = missing_data ? 1.0 : gsl_vector_get(emissionProbs[i_Ttot][1], aij);
      foreach(bkl; 0 .. msmc.nrStates) {
        auto val = gsl_matrix_get(transitionMatrix, aij, bkl) * e;
        gsl_matrix_set(ret[0], aij, bkl, val);
      }
    }

    foreach(distance; 1 .. maxDistance) {
      // if(distance % 100 == 0)
      //   stderr.writeln("distance ", distance);
      gsl_matrix_set_zero(ret[distance]);
      
      gsl_blas_dgemm_checked(CBLAS_TRANSPOSE_t.CblasNoTrans, CBLAS_TRANSPOSE_t.CblasNoTrans,
                     1.0, ret[0], ret[distance - 1], 0.0, ret[distance]);
      
    }
  }
  
  private void computeBackwardPropagators(gsl_matrix*[] ret,
                                  bool missing_data, size_t i_Ttot, size_t maxDistance) const
  {
    foreach(aij; 0 .. msmc.nrStates) {
      double e = missing_data ? 1.0 : gsl_vector_get(emissionProbs[i_Ttot][1], aij);
      foreach(bkl; 0 .. msmc.nrStates) {
        auto val = msmc.transitionRate.transitionProbability(aij, bkl) * e;
        gsl_matrix_set(ret[0], aij, bkl, val);
      }
    }
    
    foreach(distance; 1 .. maxDistance) {
      // if(distance % 100 == 0)
      //   stderr.writeln("distance ", distance);
      gsl_matrix_set_zero(ret[distance]);

      gsl_blas_dgemm_checked(CBLAS_TRANSPOSE_t.CblasNoTrans, CBLAS_TRANSPOSE_t.CblasNoTrans,
        1.0, ret[distance - 1], ret[0], 0.0, ret[distance]);
    }
  }
  
  private double fullE(in SegSite_t segsite, size_t aij) const {
    double ret = 0.0;
    foreach(o; segsite.obs) {
      ret += gsl_vector_get(emissionProbs[segsite.i_Ttot][o], aij);
    }
    ret /= cast(double)segsite.obs.length;
    return ret;
  }
  
  override void propagateSingleForward(in State_t from, State_t to, 
        in SegSite_t from_segsite, in SegSite_t to_segsite) const
  in {
    assert(to_segsite.pos == from_segsite.pos + 1);
  }
  body {
    to.setZero();
    
    foreach(aij; 0 .. msmc.nrStates) {
      auto sum = 0.0;
      foreach(bkl; 0 .. msmc.nrStates) {
        sum += from.vec[bkl] * gsl_matrix_get(transitionMatrix, aij, bkl);
      }
      to.vec[aij] = fullE(to_segsite, aij) * sum;
    }
    // gsl_blas_dgemv_checked(CBLAS_TRANSPOSE_t.CblasNoTrans, 1.0, transitionMatrix, from.getConstVector(), 0.0, to.getVector());
    // if(to_segsite.obs.length == 1) {
    //   gsl_vector_mul(to.getVector(), emissionProbs[to_segsite.i_Ttot][to_segsite.obs[0]]);
    // }
    // else {
    //   auto e_dummy = gsl_vector_alloc(msmc.nrStates);
    //   gsl_vector_memcpy(e_dummy, emissionProbs[to_segsite.i_Ttot][to_segsite.obs[0]]);
    //   foreach(o; to_segsite.obs[1..$])
    //     gsl_vector_add(e_dummy, emissionProbs[to_segsite.i_Ttot][o]);
    //   gsl_vector_scale(e_dummy, 1.0 / to_segsite.obs.length);
    //   gsl_vector_mul(to.getVector(), e_dummy);
    //   gsl_vector_free(e_dummy);
    // }
  }
  
  override void propagateSingleBackward(in State_t to, State_t from,
            in SegSite_t to_segsite, in SegSite_t from_segsite) const
  in {
    assert(to_segsite.pos == from_segsite.pos + 1);
  }
  body {
    foreach(bkl; 0 .. msmc.nrStates) {
      auto sum = 0.0;
      foreach(aij; 0 .. msmc.nrStates) {
        sum += to.vec[aij] * fullE(to_segsite, aij) * gsl_matrix_get(transitionMatrix, aij, bkl);
      }
      from.vec[bkl] = sum;
    }

    // auto e_dummy = gsl_vector_alloc(msmc.nrStates);
    // gsl_vector_memcpy(e_dummy, emissionProbs[to_segsite.i_Ttot][to_segsite.obs[0]]);
    // if(to_segsite.obs.length > 1) {
    //   foreach(o; to_segsite.obs[1..$])
    //     gsl_vector_add(e_dummy, emissionProbs[to_segsite.i_Ttot][o]);
    //   gsl_vector_scale(e_dummy, 1.0 / to_segsite.obs.length);
    // }
    // gsl_vector_mul(e_dummy, to.getConstVector());
    // gsl_blas_dgemv_checked(CBLAS_TRANSPOSE_t.CblasTrans, 1.0, transitionMatrix, e_dummy, 0.0, from.getVector());
    // gsl_vector_free(e_dummy);
  }
  
  override void propagateMultiForward(in State_t from, State_t to,
        in SegSite_t from_segsite, in SegSite_t to_segsite) const
  in {
    assert(to_segsite.pos > from_segsite.pos);
    assert(to_segsite.obs[0] < 2);
  }
  body {
    auto dist = to_segsite.pos - from_segsite.pos;
    foreach(aij; 0 .. msmc.nrStates) {
      if(to_segsite.obs[0] == 0) {
        auto prop = forwardPropagatorsMissing[dist - 1];
        auto sum = 0.0;
        foreach(bkl; 0 .. msmc.nrStates) {
          sum += from.vec[bkl] * gsl_matrix_get(prop, aij, bkl);
        }
        to.vec[aij] = sum;
      }
      else {
        auto prop = forwardPropagators[to_segsite.i_Ttot][dist - 1];
        auto sum = 0.0;
        foreach(bkl; 0 .. msmc.nrStates) {
          sum += from.vec[bkl] * gsl_matrix_get(prop, aij, bkl);
        }
        to.vec[aij] = sum;
      }
    }

  }
  
  override void propagateMultiBackward(in State_t to, State_t from,
        in SegSite_t to_segsite, in SegSite_t from_segsite) const
  in {
    assert(to_segsite.pos > from_segsite.pos);
    assert(to_segsite.obs[0] < 2);
  }
  body {  
    auto dist = to_segsite.pos - from_segsite.pos;
    foreach(bkl; 0 .. msmc.nrStates) {
      if(to_segsite.obs[0] == 0) {
        auto prop = backwardPropagatorsMissing[dist - 1];
        auto sum = 0.0;
        foreach(aij; 0 .. msmc.nrStates) {
          sum += to.vec[aij] * gsl_matrix_get(prop, aij, bkl);
        }
        from.vec[bkl] = sum;
      }
      else {
        auto prop = backwardPropagators[to_segsite.i_Ttot][dist - 1];
        auto sum = 0.0;
        foreach(aij; 0 .. msmc.nrStates) {
          sum += to.vec[aij] * gsl_matrix_get(prop, aij, bkl);
        }
        from.vec[bkl] = sum;
      }
    }

  }
  
  override string toString() const {
    return "PropagationCoreNaive";
  }
  
  override const(MSMCmodel) getMSMC() const {
    return msmc;
  }
  
  override @property size_t forwardStateSize() const {
    return msmc.nrStates;
  }

  override @property size_t backwardStateSize() const {
    return msmc.nrStates;
  }
  
  override State_t newForwardState() const {
    return new State_t(msmc.nrStates, 0, 0);
  }

  override State_t newBackwardState() const {
    return new State_t(msmc.nrStates, 0, 0);
  }

  override State_t newForwardState(StateVecAllocator stateAllocator) const {
    return new State_t(msmc.nrStates, 0, 0, stateAllocator);
  }

  override State_t newBackwardState(StateVecAllocator stateAllocator) const {
    return new State_t(msmc.nrStates, 0, 0, stateAllocator);
  }
  
  override void initialState(State_t s) const {
    foreach(aij; 0 .. msmc.nrStates) {
      auto val = msmc.transitionRate.equilibriumProbability(aij);
      s.vec[aij] = val;
    }
  }
  
  override void setState(State_t s, double x, in SegSite_t segsite) const {
    foreach(aij; 0 .. msmc.nrStates)
      s.vec[aij] = x;
  }
  
  override void getTransitionExpectation(State_t f, State_t b,
      in SegSite_t to_segsite, double[][] ret) const
  {
    foreach(au; 0 .. msmc.marginalIndex.nrMarginals)
      ret[au][] = 0.0;
    foreach(aij; 0 .. msmc.nrStates) {
      auto au = msmc.marginalIndex.getMarginalIndexFromIndex(aij);

      auto e = 0.0;
      foreach(o; to_segsite.obs)
        e += gsl_vector_get(emissionProbs[to_segsite.i_Ttot][o], aij);
      e /= cast(double)to_segsite.obs.length;
      foreach(bkl; 0 .. msmc.nrStates) {
        auto bv = msmc.marginalIndex.getMarginalIndexFromIndex(bkl);
        ret[au][bv] += f.vec[bkl] * gsl_matrix_get(transitionMatrix, aij, bkl) *
                       b.vec[aij] * e;
      }
    }
  }
  
  override @property size_t maxDistance() const {
    return cast(size_t)forwardPropagators[0].length;
  }
  
}

