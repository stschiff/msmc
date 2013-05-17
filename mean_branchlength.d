unittest {
  assert(fact(2) == 2);
  assert(fact(3) == 6);
  assert(fact(4) == 24);
  assert(fact(5) == 120);
}

size_t fact(size_t n)
{
  if(n <= 1) return 1;
  auto ret = 1;
  foreach(i; 2 .. n + 1) {
    ret *= i;
  }
  return ret;
}

unittest {
  assert(isNaN(gammInt(0)));
  assert(gammInt(1) == 1);
  assert(gammInt(2) == 1);
  assert(gammInt(3) == 2);
  assert(gammInt(4) == 6);
}

double gammInt(size_t n) {
  if(n == 0)
    return double.nan;
  return cast(double)fact(n - 1);
}

private double cCoeff(int M, int m, int j) const {
  double ret = (-1.0) ^^ cast(double)(j - m);
  ret *= (2 * j - 1) / cast(double)(fact(m) * fact(j - m));
  ret *= gammInt(m + j - 1) / cast(double)gammInt(m);
  ret *= gammInt(M) / cast(double)gammInt(M + j);
  ret *= gammInt(M + 1) / cast(double)gammInt(M + 1 - j);
  return ret;
}
  
double mean_Ttot(int M, double t, int a) const {
  auto Ttot = 0.0;
  foreach(m; 2 .. M + 1) {
    foreach(j; m .. M + 1) {
      auto firstTerm = (1.0 - exp(-mOver2(j) * avgLambdaVec[a] * (modelSpecs.timeBoundaries[a + 1] - t))) / (mOver2(j) * avgLambdaVec[a]);
      auto secondTerm = 0.0;
      foreach(g; a + 1 .. T) {
        secondTerm += LIntegral(avgLambdaVec, t, a, modelSpecs.timeBoundaries[g], g) ^^ mOver2(j) * (1.0 - exp(-mOver2(j) * avgLambdaVec[g] * delta[g])) / (mOver2(j) * avgLambdaVec[g]);
      }
      auto c = cCoeff(M, m, j);
      Ttot += m * c * (firstTerm + secondTerm);
    }
  }
  return Ttot;
}
