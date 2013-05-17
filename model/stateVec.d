import std.conv;
import std.algorithm;
import gsl_matrix_vector;
import stateVecAllocator;

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
