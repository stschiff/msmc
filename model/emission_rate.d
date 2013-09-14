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
 
module model.emission_rate;
import std.stdio;
import std.math;
import std.exception;
import std.conv;
import std.algorithm;
import std.range;
import model.triple_index;

class EmissionRate {
  size_t nrHaplotypes;
  double mu;
  bool directedEmissions;
  
  this(size_t nrHaplotypes, double mu, bool directedEmissions) {
    enforce(mu > 0, "need positive mutation rate");
    enforce(nrHaplotypes >= 2, "need at least two haplotypes");
    this.mu = mu;
    this.nrHaplotypes = nrHaplotypes;
    this.directedEmissions = directedEmissions;
  }
  
  int getEmissionId(string alleles, size_t ind1, size_t ind2) const {
    if(!directedEmissions) {
      auto count_0 = count(alleles, '0');
      auto count_1 = nrHaplotypes - count_0;
      if(count_0 == 0 || count_1 == 0)
        return 0;
      if(count_0 == 1 || count_1 == 1) {
        if(alleles[ind1] != alleles[ind2])
          return 1;
        else
          return cast(int)nrHaplotypes - 1;
      }
      else {
        if(alleles[ind1] != alleles[ind2])
          return -1;
        else
          return cast(int)count(alleles, alleles[ind1]);
      }
    }
    else {
      return getEmissionIdDirected(alleles, ind1, ind2);
    }
  }
  
  int getEmissionIdDirected(string alleles, size_t ind1, size_t ind2) const {
    auto count_0 = count(alleles, '0');
    auto count_1 = nrHaplotypes - count_0;
    if(count_1 == 0)
      return 0;
    if(count_1 == 1) {
      if(alleles[ind1] != alleles[ind2])
        return 1;
      else
        return cast(int)nrHaplotypes;
    }
    if(alleles[ind1] != alleles[ind2])
      return -1;
    if(alleles[ind1] == '1')
      return cast(int)count_1;
    else
      return cast(int)(nrHaplotypes + count_1 - 1);
  }
  
  size_t getNrEmissionIds() const {
    if(directedEmissions)
      return 2 * nrHaplotypes - 2;
    else
      return nrHaplotypes;
  }
  
  double emissionProb(int emissionId, double t, double tLeaf) const {
    if(directedEmissions)
      return emissionProbDirected(emissionId, t, tLeaf);
    if(nrHaplotypes == 2) {
      if(emissionId == 0)
        return exp(-2.0 * mu * t);
      else
        return 1.0 - exp(-2.0 * mu * t);
    }
    
    auto tLower = t * nrHaplotypes;
    if(tLeaf < tLower)
      tLeaf = tLower;

    auto tLeafUpper = tLeaf - tLower;
    auto tUpper = tLeafUpper + mutationTreeLength(nrHaplotypes - 1, 1) +
                  mutationTreeLength(nrHaplotypes - 1, nrHaplotypes - 2);
    if(nrHaplotypes > 4)
      tUpper += 2.0 * iota(2, nrHaplotypes - 2).map!"1.0/a"().reduce!"a+b"();
    auto tTot = tLower + tUpper;
    
    if(emissionId < 0)
      return 0.0;
    if(emissionId == 0) {
      return exp(-mu * tTot);
    }
    if(emissionId == 1)
      return (1.0 - exp(-mu * tTot)) * t / tTot;
    if(emissionId == nrHaplotypes - 1)
      return (1.0 - exp(-mu * tTot)) * (tLeaf - 2.0 * t) / tTot / (nrHaplotypes - 2.0);
    else {
      return (1.0 - exp(-mu * tTot)) * (mutationTreeLength(nrHaplotypes - 1, emissionId - 1) + 
                                        mutationTreeLength(nrHaplotypes - 1, emissionId)) / tTot;
      // auto num = 1.0 + 1.0;
      // if(nrHaplotypes > 4)
      //   foreach(k; 2 .. nrHaplotypes - 2)
      //     num += binomial(nrHaplotypes - 1, k);
      // return (1.0 - exp(-mu * tTot)) * 2.0 * (tTot - tLeaf) / tTot / num;
    }
  }

  double emissionProbDirected(int emissionId, double t, double tLeaf) const {
    if(nrHaplotypes == 2) {
      if(emissionId == 0)
        return exp(-2.0 * mu * t);
      else
        return 1.0 - exp(-2.0 * mu * t);
    }
    
    auto tLower = t * nrHaplotypes;
    if(tLeaf < tLower)
      tLeaf = tLower;

    auto tLeafUpper = tLeaf - tLower;
    auto tUpper = tLeafUpper + mutationTreeLength(nrHaplotypes - 1, 1);
    tUpper += 2.0 * iota(2, nrHaplotypes - 1).map!"1.0/a"().reduce!"a+b"();
    auto tTot = tLower + tUpper;
    
    if(emissionId < 0)
      return 0.0;
    if(emissionId == 0) {
      return exp(-mu * tTot);
    }
    if(emissionId == 1)
      return (1.0 - exp(-mu * tTot)) * t / tTot;
    if(emissionId == nrHaplotypes)
      return (1.0 - exp(-mu * tTot)) * (tLeaf - 2.0 * t) / tTot / (nrHaplotypes - 2.0);
    else {
      if(emissionId < nrHaplotypes)
        return (1.0 - exp(-mu * tTot)) * mutationTreeLength(nrHaplotypes - 1, emissionId - 1) / tTot;
      else {
        auto freq = emissionId - nrHaplotypes + 1;
        return (1.0 - exp(-mu * tTot)) * mutationTreeLength(nrHaplotypes - 1, freq) / tTot;
      }
    }
  }
}

double mutationTreeLength(size_t m, size_t freq)
in {
  assert(m >= freq);
  assert(m > 0);
  assert(freq > 0);
}
out(res) {
  assert(res > 0.0);
  assert(res < 1.0);
}
body {
  return 2.0 / freq * (1.0 / binomial(m, freq));
}


double binomial(size_t m, size_t k) {
  auto prod = 1.0;
  foreach(i; 1 .. k + 1) {
    prod *= cast(double)(m - (k - i)) / i;
  }
  return prod;
}

unittest {
  assert(binomial(6,2) == 15);
  assert(binomial(4,2) == 6);
  assert(binomial(6,5) == 6);
  assert(binomial(6,6) == 1);
  assert(binomial(4,4) == 1);
}

unittest {
  writeln("test emissionRate.emissionType");
  auto emissionRate = new EmissionRate(4, 0.01, false);
  assert(emissionRate.getEmissionId("0000", 0, 1) == 0);
  assert(emissionRate.getEmissionId("0001", 0, 1) == 3);
  assert(emissionRate.getEmissionId("0011", 0, 1) == 2);
  assert(emissionRate.getEmissionId("0100", 0, 1) == 1);
  assert(emissionRate.getEmissionId("1000", 0, 1) == 1);
  assert(emissionRate.getEmissionId("1001", 0, 1) == -1);
}

unittest {
  writeln("test emissionRate.emissionProb");
  auto emissionRate = new EmissionRate(4, 0.01, false);
  assert(emissionRate.emissionProb(-1, 1.0, 2.0) == 0.0);
  assert(emissionRate.emissionProb(-1, 0.0, 2.0) == 0.0);
  assert(emissionRate.emissionProb(-1, 0.5, 2.0) == 0.0);
  assert(emissionRate.emissionProb(-1, 1.5, 2.0) == 0.0);
  assert(emissionRate.emissionProb(0, 1.0, 2.0) > 0.0);
  assert(emissionRate.emissionProb(0, 1.0, 2.0) < 1.0);
  assert(emissionRate.emissionProb(0, 0.0, 2.0) > 0.0);
  assert(emissionRate.emissionProb(0, 0.0, 2.0) < 1.0);
  assert(emissionRate.emissionProb(0, 0.5, 2.0) > 0.0);
  assert(emissionRate.emissionProb(0, 0.5, 2.0) < 1.0);
  assert(emissionRate.emissionProb(0, 1.5, 2.0) > 0.0);
  assert(emissionRate.emissionProb(0, 1.5, 2.0) < 1.0);
  assert(emissionRate.emissionProb(2, 0.0, 2.0) > 0.0);
  assert(emissionRate.emissionProb(2, 0.0, 2.0) < 1.0);
  assert(emissionRate.emissionProb(2, 0.5, 2.0) > 0.0);
  assert(emissionRate.emissionProb(2, 0.5, 2.0) < 1.0);
  assert(emissionRate.emissionProb(1, 0.5, 2.0) > 0.0);
  assert(emissionRate.emissionProb(1, 0.5, 2.0) < 1.0);
}