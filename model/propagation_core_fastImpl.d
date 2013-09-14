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

module model.propagation_core_fastImpl;
import std.stdio;
import std.algorithm;
import std.string;
import std.conv;
import std.math;
import std.exception;
import model.msmc_hmm;
import model.data;
import model.gsl_matrix_vector;
import model.propagation_core;
import model.msmc_model;
import model.stateVec;
import model.stateVecAllocator;

double gsl_matrix_get_checked(const gsl_matrix* m, size_t i, size_t j) {
  assert(i < (*m).size1);
  assert(j < (*m).size2);
  return gsl_matrix_get(m, i, j);
}
  
double gsl_vector_get_checked(const gsl_vector* v, size_t i) {
  assert(i < (*v).size);
  return gsl_vector_get(v, i);
}
  
class PropagationCoreFast : PropagationCore {
  gsl_vector*[][] propagatorsDiag;
  gsl_matrix*[][] forwardPropagatorsOffDiag, backwardPropagatorsOffDiag;
  gsl_vector*[] propagatorsDiagMissing;
  gsl_matrix*[] forwardPropagatorsOffDiagMissing, backwardPropagatorsOffDiagMissing;
  
  double[][][] emissionProbs; // first index: i_ttot, second index: obs, third: state
  double[][] emissionProbsMarginal;
  double[] transitionMatrixQ1;
  double[][] transitionMatrixQ2;

  const MSMCmodel msmc;
  
  this(in MSMCmodel msmc, size_t maxDistance) {
    this.msmc = msmc;
    enforce(maxDistance > 0);

    auto allele_order = canonicalAlleleOrder(msmc.nrHaplotypes);
    
    emissionProbs = new double[][][](msmc.nrTtotIntervals, allele_order.length + 1, msmc.nrStates);
    foreach(tt; 0 .. msmc.nrTtotIntervals) {
      foreach(i; 0 .. allele_order.length + 1) {
        foreach(aij; 0 .. msmc.nrStates) {
          if(i == 0)
            emissionProbs[tt][i][aij] = 1.0; // missing data
          else
            emissionProbs[tt][i][aij] = msmc.emissionProb(allele_order[i - 1], aij, tt);
        }
      }
    }

    emissionProbsMarginal = new double[][](msmc.nrTtotIntervals, msmc.nrMarginals);
    foreach(tt; 0 .. msmc.nrTtotIntervals) {
      foreach(au; 0 .. msmc.nrMarginals) {
        auto index = msmc.marginalIndex.getIndexFromMarginalIndex(au);
        auto time = msmc.marginalIndex.getTripleFromIndex(index).time;
        emissionProbsMarginal[tt][au] = msmc.emissionProbHom(time, tt);
      }
    }
    
    transitionMatrixQ1 = new double[msmc.nrMarginals];
    transitionMatrixQ2 = new double[][](msmc.nrMarginals, msmc.nrMarginals);
    foreach(au; 0 .. msmc.nrMarginals) {
      transitionMatrixQ1[au] = msmc.transitionRate.transitionProbabilityQ1(au);
      foreach(bv; 0 .. msmc.nrMarginals) {
        transitionMatrixQ2[au][bv] = msmc.transitionRate.transitionProbabilityQ2(au, bv);
      }
    }
    // stderr.writeln("pre-computing propagators for homozygous segments");
    propagatorsDiag = new gsl_vector*[][](msmc.nrTtotIntervals, maxDistance);
    forwardPropagatorsOffDiag = new gsl_matrix*[][](msmc.nrTtotIntervals, maxDistance);
    backwardPropagatorsOffDiag = new gsl_matrix*[][](msmc.nrTtotIntervals, maxDistance);
      
    foreach(i; 0 .. msmc.nrTtotIntervals) {
      foreach(d; 0 .. maxDistance) {
        propagatorsDiag[i][d] = gsl_vector_alloc(msmc.nrMarginals);
        forwardPropagatorsOffDiag[i][d] = gsl_matrix_alloc(msmc.nrMarginals, msmc.nrMarginals);
        backwardPropagatorsOffDiag[i][d] = gsl_matrix_alloc(msmc.nrMarginals, msmc.nrMarginals);
      }
      // stderr.writeln("pre-computing propagators for homozygous segments, i_Ttot=", i);
      computePropagatorsDiag(propagatorsDiag[i], false, i, maxDistance);
      computeForwardPropagatorsOffDiag(forwardPropagatorsOffDiag[i], false, i, maxDistance);
      computeBackwardPropagatorsOffDiag(backwardPropagatorsOffDiag[i], false, i, maxDistance);
    }
    propagatorsDiagMissing = new gsl_vector*[maxDistance];
    forwardPropagatorsOffDiagMissing = new gsl_matrix*[maxDistance];
    backwardPropagatorsOffDiagMissing = new gsl_matrix*[maxDistance];
    foreach(d; 0 .. maxDistance) {
      propagatorsDiagMissing[d] = gsl_vector_alloc(msmc.nrMarginals);
      forwardPropagatorsOffDiagMissing[d] = gsl_matrix_alloc(msmc.nrMarginals, msmc.nrMarginals);
      backwardPropagatorsOffDiagMissing[d] = gsl_matrix_alloc(msmc.nrMarginals, msmc.nrMarginals);
    }
    // stderr.writeln("pre-computing propatagors for missing data segments");
    computePropagatorsDiag(propagatorsDiagMissing, true, 0, maxDistance);
    computeForwardPropagatorsOffDiag(forwardPropagatorsOffDiagMissing, true, 0, maxDistance);
    computeBackwardPropagatorsOffDiag(backwardPropagatorsOffDiagMissing, true, 0, maxDistance);
    
  }
  
  ~this() {
    foreach(i; 0 .. propagatorsDiag.length) {
      foreach(d; 0 .. propagatorsDiag[i].length) {
        gsl_vector_free(propagatorsDiag[i][d]);
        gsl_matrix_free(forwardPropagatorsOffDiag[i][d]);
        gsl_matrix_free(backwardPropagatorsOffDiag[i][d]);
      }
    }
    foreach(d; 0 .. propagatorsDiagMissing.length) {
      gsl_vector_free(propagatorsDiagMissing[d]);
      gsl_matrix_free(forwardPropagatorsOffDiagMissing[d]);
      gsl_matrix_free(backwardPropagatorsOffDiagMissing[d]);
    }
  } 
  
  private void computePropagatorsDiag(gsl_vector*[] ret, bool missing_data, size_t i, size_t maxDistance) const {
    foreach(au; 0 .. msmc.nrMarginals) {
      auto e = missing_data ? 1.0 : emissionProbsMarginal[i][au];
      gsl_vector_set(ret[0], au, transitionMatrixQ1[au] * e);
    }
    
    foreach(distance; 1 .. maxDistance) {
      gsl_vector_memcpy(ret[distance], ret[0]);
      gsl_vector_mul(ret[distance], ret[distance - 1]);
    }
  }
  
  private void computeForwardPropagatorsOffDiag(
      gsl_matrix*[] ret, bool missing_data, size_t i, size_t maxDistance) const
  {
    auto q2cardinal = gsl_matrix_alloc(msmc.nrMarginals, msmc.nrMarginals);
    auto e = new double[msmc.nrMarginals];
    
    foreach(au; 0 .. msmc.nrMarginals) {
      e[au] = missing_data ? 1.0 : emissionProbsMarginal[i][au];
      foreach(bv; 0 .. msmc.nrMarginals) {
        gsl_matrix_set(q2cardinal, au, bv, transitionMatrixQ2[au][bv] *
            cast(double)msmc.marginalIndex.getDegeneracyForMarginalIndex(bv));
        gsl_matrix_set(ret[0], au, bv, transitionMatrixQ2[au][bv] * e[au]);
      }
    }
    
    foreach(distance; 1 .. maxDistance) {
      gsl_blas_dgemm_checked(CBLAS_TRANSPOSE_t.CblasNoTrans, CBLAS_TRANSPOSE_t.CblasNoTrans,
                     1.0, q2cardinal, ret[distance - 1], 0.0, ret[distance]);
      
      foreach(bv; 0 .. msmc.nrMarginals) {
        double g1;
        // distance here means distance - 1, since we use 0-based indices
        g1 = (transitionMatrixQ1[bv] * e[bv]) ^^ distance;
        foreach(au; 0 .. msmc.nrMarginals) {
          auto val = gsl_matrix_get(ret[distance - 1], au, bv) * transitionMatrixQ1[au] +
                     g1 * transitionMatrixQ2[au][bv] + gsl_matrix_get(ret[distance], au, bv);
          gsl_matrix_set(ret[distance], au, bv, val * e[au]);
        }
      }
    }
    gsl_matrix_free(q2cardinal);
  }
  
  private void computeBackwardPropagatorsOffDiag(
      gsl_matrix*[] ret, bool missing_data, size_t i, size_t maxDistance) const
  {

    auto q2cardinalEmission = gsl_matrix_alloc(msmc.nrMarginals, msmc.nrMarginals);
    auto e = new double[msmc.nrMarginals];
    
    foreach(au; 0 .. msmc.nrMarginals) {
      auto index = msmc.marginalIndex.getIndexFromMarginalIndex(au);
      auto triple = msmc.marginalIndex.getTripleFromIndex(index);
      e[au] = missing_data ? 1.0 : msmc.emissionProbHom(triple.time, i);
    }
    
    foreach(au; 0 .. msmc.nrMarginals) {
      foreach(bv; 0 .. msmc.nrMarginals) {
        auto val = e[au] * transitionMatrixQ2[au][bv] * cast(double)msmc.marginalIndex.getDegeneracyForMarginalIndex(au);
        gsl_matrix_set(q2cardinalEmission, au, bv, val);
        gsl_matrix_set(ret[0], au, bv, transitionMatrixQ2[au][bv] * e[au]);
      }
    }
    
    foreach(distance; 1 .. maxDistance) {
      gsl_blas_dgemm_checked(CBLAS_TRANSPOSE_t.CblasNoTrans, CBLAS_TRANSPOSE_t.CblasNoTrans,
                     1.0, ret[distance - 1], q2cardinalEmission, 0.0, ret[distance]);
      
      foreach(au; 0 .. msmc.nrMarginals) {
        double h1;
        // distance here means distance - 1, since we use 0-based indices
        h1 = (transitionMatrixQ1[au] * e[au]) ^^ distance;
        foreach(bv; 0 .. msmc.nrMarginals) {
          auto val = gsl_matrix_get(ret[distance - 1], au, bv) * transitionMatrixQ1[bv] * e[bv] +
              h1 * transitionMatrixQ2[au][bv] * e[au] + gsl_matrix_get(ret[distance], au, bv);
          gsl_matrix_set(ret[distance], au, bv, val);
        }
      }
    }
    gsl_matrix_free(q2cardinalEmission);
  }
  
  private double fullE(in SegSite_t segsite, size_t aij) const {
    double ret = 0.0;
    foreach(o; segsite.obs) {
      ret += emissionProbs[segsite.i_Ttot][o][aij];
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
      auto au = msmc.marginalIndex.getMarginalIndexFromIndex(aij);
      auto sum = 0.0;
      foreach(bv; 0 .. msmc.nrMarginals) {
        sum += from.vecMarginal[bv] * transitionMatrixQ2[au][bv];
      }
      to.vec[aij] = fullE(to_segsite, aij) * (sum + from.vec[aij] * transitionMatrixQ1[au]);
      to.vecMarginal[au] = to.vecMarginal[au] + to.vec[aij];
    }
  }

  override void propagateSingleBackward(in State_t to, State_t from,
            in SegSite_t to_segsite, in SegSite_t from_segsite) const
  in {
    assert(to_segsite.pos == from_segsite.pos + 1);
  }
  body {
    from.setZero();
    
    foreach(bkl; 0 .. msmc.nrStates) {
      auto bv = msmc.marginalIndex.getMarginalIndexFromIndex(bkl);
      auto sum = 0.0;
      foreach(au; 0 .. msmc.nrMarginals) {
        sum += to.vecMarginalEmission[au] * transitionMatrixQ2[au][bv];
      }
      from.vec[bkl] = sum + to.vec[bkl] * fullE(to_segsite, bkl) * transitionMatrixQ1[bv];
      from.vecMarginal[bv] = from.vecMarginal[bv] + from.vec[bkl];
      from.vecMarginalEmission[bv] = from.vecMarginalEmission[bv] + from.vec[bkl] * fullE(from_segsite, bkl);
    }
  }
  
  override void propagateMultiForward(in State_t from, State_t to,
        in SegSite_t from_segsite, in SegSite_t to_segsite) const
  in {
    assert(to_segsite.pos > from_segsite.pos);
    assert(to_segsite.obs[0] < 2);
  }
  body {
    auto dist = to_segsite.pos - from_segsite.pos;
    to.setZero();
    foreach(aij; 0 .. msmc.nrStates) {
      auto au = msmc.marginalIndex.getMarginalIndexFromIndex(aij);
      if(to_segsite.obs[0] == 0) {
        auto prop = forwardPropagatorsOffDiagMissing[dist - 1];
        auto propDiag = propagatorsDiagMissing[dist - 1];
        auto sum = 0.0;
        foreach(bv; 0 .. msmc.nrMarginals) {
          sum += from.vecMarginal[bv] * gsl_matrix_get_checked(prop, au, bv);
        }
        to.vec[aij] = sum + from.vec[aij] * gsl_vector_get_checked(propDiag, au);
      }
      else {
        auto prop = forwardPropagatorsOffDiag[to_segsite.i_Ttot][dist - 1];
        auto propDiag = propagatorsDiag[to_segsite.i_Ttot][dist - 1];
        auto sum = 0.0;
        foreach(bv; 0 .. msmc.nrMarginals) {
          sum += from.vecMarginal[bv] * gsl_matrix_get_checked(prop, au, bv);
        }
        to.vec[aij] = sum + from.vec[aij] * gsl_vector_get_checked(propDiag, au);
      }
      to.vecMarginal[au] = to.vecMarginal[au] + to.vec[aij];
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
    from.setZero();
    foreach(bkl; 0 .. msmc.nrStates) {
      auto bv = msmc.marginalIndex.getMarginalIndexFromIndex(bkl);
      if(to_segsite.obs[0] == 0) {
        auto prop = backwardPropagatorsOffDiagMissing[dist - 1];
        auto propDiag = propagatorsDiagMissing[dist - 1];
        auto sum = 0.0;
        foreach(au; 0 .. msmc.nrMarginals) {
          sum += to.vecMarginal[au] * gsl_matrix_get_checked(prop, au, bv);
        }
        from.vec[bkl] = sum + to.vec[bkl] * gsl_vector_get_checked(propDiag, bv);
      }
      else {
        auto prop = backwardPropagatorsOffDiag[to_segsite.i_Ttot][dist - 1];
        auto propDiag = propagatorsDiag[to_segsite.i_Ttot][dist - 1];
        auto sum = 0.0;
        foreach(au; 0 .. msmc.nrMarginals) {
          sum += to.vecMarginal[au] * gsl_matrix_get_checked(prop, au, bv);
        }
        from.vec[bkl] = sum + to.vec[bkl] * gsl_vector_get_checked(propDiag, bv);
      }

      from.vecMarginal[bv] = from.vecMarginal[bv] + from.vec[bkl];
      from.vecMarginalEmission[bv] = from.vecMarginalEmission[bv] + from.vec[bkl] * fullE(from_segsite, bkl);
    }
  }
  
  override string toString() const {
    return "PropagationCoreFast";
  }
  
  override const(MSMCmodel) getMSMC() const {
    return msmc;
  }
  
  override @property size_t forwardStateSize() const {
    return msmc.nrStates + msmc.nrMarginals;
  }

  override @property size_t backwardStateSize() const {
    return msmc.nrStates + 2 * msmc.nrMarginals;
  }
  
  override State_t newForwardState() const {
    return new State_t(msmc.nrStates, msmc.nrMarginals, 0);
  }

  override State_t newBackwardState() const {
    return new State_t(msmc.nrStates, msmc.nrMarginals, msmc.nrMarginals);
  }

  override State_t newForwardState(StateVecAllocator stateAllocator) const {
    return new State_t(msmc.nrStates, msmc.nrMarginals, 0, stateAllocator);
  }

  override State_t newBackwardState(StateVecAllocator stateAllocator) const {
    return new State_t(msmc.nrStates, msmc.nrMarginals, msmc.nrMarginals, stateAllocator);
  }

  override void initialState(State_t s) const {
    s.setZero();
    foreach(aij; 0 .. msmc.nrStates) {
      s.vec[aij] = msmc.transitionRate.equilibriumProbability(aij);
      auto au = msmc.marginalIndex.getMarginalIndexFromIndex(aij);
      s.vecMarginal[au] += s.vec[aij];
    }
  }
  
  override void setState(State_t s, double x, in SegSite_t segsite) const {
    s.setZero();
    foreach(aij; 0 .. msmc.nrStates)
      s.vec[aij] = x;
    
    foreach(aij; 0 .. msmc.nrStates) {
      auto au = msmc.marginalIndex.getMarginalIndexFromIndex(aij);
      s.vecMarginal[au] += x;
      if(s.vecMarginalEmission.length > 0)
        s.vecMarginalEmission[au] += x * fullE(segsite, aij);
    }
  }
  
  override void getTransitionExpectation(State_t f, State_t b,
      in SegSite_t to_segsite, double[] eVec, double[][] eMat) const
  {
    foreach(bv; 0 .. msmc.nrMarginals) {
      foreach(au; 0 .. msmc.nrMarginals) {
        eMat[au][bv] =
            f.vecMarginal[bv] * transitionMatrixQ2[au][bv] * b.vecMarginalEmission[au];
      }
    }
    eVec[] = 0.0;
    foreach(aij; 0 .. msmc.nrStates) {
      auto au = msmc.marginalIndex.getMarginalIndexFromIndex(aij);
      eMat[au][au] -= f.vec[aij] * transitionMatrixQ2[au][au] * fullE(to_segsite, aij) * b.vec[aij];
      if(eMat[au][au] < 0.0) { // this just corrects tiny numerical errors from the addition-subtraction here.
        assert(-eMat[au][au] < 1.0e-20);
        eMat[au][au] = 0.0;
      }
      eVec[au] += f.vec[aij] * (transitionMatrixQ1[au] + transitionMatrixQ2[au][au]) * fullE(to_segsite, aij) * 
                  b.vec[aij];
    }
  }
  
  override @property size_t maxDistance() const {
    return cast(size_t)forwardPropagatorsOffDiag[0].length;
  }
}

unittest {
  writeln("testing propagationCoreFast and propagationCoreNaive propagators");

  import model.propagation_core_naiveImpl;
  auto lambdaVec = new double[30];
  lambdaVec[] = 1.0;
  auto msmc = new MSMCmodel(0.01, 0.001, [0U, 0, 1, 1], lambdaVec, 10, 4, false);
  auto lvl = 1.0e-8;
  
  auto maxDist = 10U;

  auto propagationCoreNaive = new PropagationCoreNaive(msmc, maxDist);
  auto propagationCoreFast = new PropagationCoreFast(msmc, maxDist);

  auto d = 1U;
  auto n_tTot = 0;
  foreach(aij; 0 .. msmc.nrStates) {
    auto au = msmc.marginalIndex.getMarginalIndexFromIndex(aij);
    foreach(bkl; 0 .. msmc.nrStates) {
      auto bv = msmc.marginalIndex.getMarginalIndexFromIndex(bkl);
      auto gConstr = gsl_matrix_get(propagationCoreFast.forwardPropagatorsOffDiag[n_tTot][d], au, bv);
      auto hConstr = gsl_matrix_get(propagationCoreFast.backwardPropagatorsOffDiag[n_tTot][d], au, bv);
      auto gConstrM = gsl_matrix_get(propagationCoreFast.forwardPropagatorsOffDiagMissing[d], au, bv);
      auto hConstrM = gsl_matrix_get(propagationCoreFast.backwardPropagatorsOffDiagMissing[d], au, bv);
      if(aij == bkl) {
        gConstr += gsl_vector_get(propagationCoreFast.propagatorsDiag[n_tTot][d], au);
        hConstr += gsl_vector_get(propagationCoreFast.propagatorsDiag[n_tTot][d], au);
        gConstrM += gsl_vector_get(propagationCoreFast.propagatorsDiagMissing[d], au);
        hConstrM += gsl_vector_get(propagationCoreFast.propagatorsDiagMissing[d], au);
      }
      auto val = gsl_matrix_get(propagationCoreNaive.forwardPropagators[n_tTot][d], aij, bkl);
      assert(approxEqual(gConstr, val, lvl, 0.0), text([gConstr, val]));
      val = gsl_matrix_get(propagationCoreNaive.backwardPropagators[n_tTot][d], aij, bkl);
      assert(approxEqual(hConstr, val, lvl, 0.0));
      val = gsl_matrix_get(propagationCoreNaive.forwardPropagatorsMissing[d], aij, bkl);
      assert(approxEqual(gConstrM,  val, lvl, 0.0));
      val = gsl_matrix_get(propagationCoreNaive.forwardPropagatorsMissing[d], aij, bkl);
      assert(approxEqual(hConstrM, val, lvl, 0.0));
    }
  }  
}

unittest {
  writeln("testing propagationCoreFast and propagationCoreNaive propagateForward ");
  import model.propagation_core_naiveImpl;
  auto lambdaVec = new double[30];
  lambdaVec[] = 1.0;
  auto msmc = new MSMCmodel(0.01, 0.001, [0U, 0, 1, 1], lambdaVec, 10, 4, false);
  auto lvl = 1.0e-8;
  auto propagationCoreNaive = new PropagationCoreNaive(msmc, 10);
  auto propagationCoreFast = new PropagationCoreFast(msmc, 10);
  
  double[] testForwardPropagation(impl_t)(impl_t impl) {
    auto dist = 10;
    auto f = impl.newForwardState();
    auto fNext = impl.newForwardState();
    auto dummy_site = new SegSite_t(1, 1, 2);
    impl.setState(f, 1.0, dummy_site);
    foreach(i; 0 .. dist) {
      auto left_site = new SegSite_t(1 + i, 1, 2);
      auto right_site = new SegSite_t(1 + i + 1, 1, 2);
      impl.propagateSingleForward(f, fNext, left_site, right_site);
      fNext.copy_into(f);
    }
    auto fSingles = impl.newForwardState();
    f.copy_into(fSingles);
    impl.setState(f, 1.0, dummy_site);
    auto left_site = new SegSite_t(1, 1, 2);
    auto right_site = new SegSite_t(1 + dist, 1, 2);
    impl.propagateMultiForward(f, fNext, left_site, right_site);
    auto fSinglesA = fSingles.vec;
    auto fNextA = fNext.vec;
    foreach(aij; 0 .. msmc.nrStates) {
      assert(
          approxEqual(fNextA[aij], fSinglesA[aij], lvl, 0.0),
          text(fNext, " ", fSingles)
      );
    }
    return fNextA;
  }  
  
  auto forwardNaive = testForwardPropagation(propagationCoreNaive);
  auto forwardFast = testForwardPropagation(propagationCoreFast);
  foreach(aij; 0 .. msmc.nrStates) {
    assert(approxEqual(forwardNaive[aij], forwardFast[aij], lvl, 0.0), text(forwardNaive, " ", forwardFast));
  }
}
  

unittest {
  writeln("testing propagationCoreFast and propagationCoreNaive propagateBackward ");
  import model.propagation_core_naiveImpl;
  auto lambdaVec = new double[30];
  lambdaVec[] = 1.0;
  auto msmc = new MSMCmodel(0.01, 0.001, [0U, 0, 1, 1], lambdaVec, 10, 4, false);
  auto lvl = 1.0e-8;
  auto propagationCoreNaive = new PropagationCoreNaive(msmc, 10);
  auto propagationCoreFast = new PropagationCoreFast(msmc, 10);
  
  double[] testBackwardPropagation(impl_t)(impl_t impl) {
    auto dist = 10;
    auto b = impl.newBackwardState();
    auto bNext = impl.newBackwardState();
    auto dummy_site = new SegSite_t(dist + 1, 1, 2);
    impl.setState(b, 1.0, dummy_site);
    foreach(i; 0 .. dist) {
      auto left_site = new SegSite_t(dist - i, 1, 2);
      auto right_site = new SegSite_t(dist + 1 - i, 1, 2);
      impl.propagateSingleBackward(b, bNext, right_site, left_site);
      bNext.copy_into(b);
    }
    auto bSingles = impl.newBackwardState();
    b.copy_into(bSingles);
    impl.setState(b, 1.0, dummy_site);
    auto left_site = new SegSite_t(1, 1, 2);
    auto right_site = new SegSite_t(dist + 1, 1, 2);
    impl.propagateMultiBackward(b, bNext, right_site, left_site);
    auto bSinglesA = bSingles.vec;
    auto bNextA = bNext.vec;
    foreach(aij; 0 .. msmc.nrStates) {
      assert(
          approxEqual(bNextA[aij], bSinglesA[aij], lvl, 0.0),
          text(bNext, ", ", bSingles)
      );
    }
    return bNextA;
  }

  auto backwardNaive = testBackwardPropagation(propagationCoreNaive);
  auto backwardFast = testBackwardPropagation(propagationCoreFast);
  foreach(aij; 0 .. msmc.nrStates) {
    assert(approxEqual(backwardNaive[aij], backwardFast[aij], lvl, 0.0), text(backwardNaive, " ", backwardFast));
  }
}
  
