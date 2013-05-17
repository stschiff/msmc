import std.conv;
import data;
import stateVec;
import stateVecAllocator;
import msmc_model;

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