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