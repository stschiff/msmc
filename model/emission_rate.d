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
 
import std.stdio;
import std.math;
import std.exception;
import std.conv;
import triple_index;

class EmissionRate {
  size_t nrHaplotypes;
  double mu;
  
  this(size_t nrHaplotypes, double mu) {
    enforce(mu > 0, "need positive mutation rate");
    enforce(nrHaplotypes >= 2, "need at least two haplotypes");
    this.mu = mu;
    this.nrHaplotypes = nrHaplotypes;
  }
  
  enum Observation_t {NoMut, MutElsewhere, MutInPair, DoubleMut}
  
  Observation_t emissionType(string alleles, size_t ind1, size_t ind2) const
    in {
      assert(ind2 > ind1);
    }
  body {
    if(alleles[ind1] == alleles[ind2]) {
      return emissionTypeHomInPair(alleles, ind1, ind2);
    }
    else { // we have a mutation in the pair
      return emissionTypeHetInPair(alleles, ind1, ind2);
    }
  }
  
  private Observation_t emissionTypeHomInPair(string alleles, size_t ind1, size_t ind2) const {
    foreach(k, a; alleles) {
      if(k != ind1 && k != ind2 && a != alleles[ind1])
        return Observation_t.MutElsewhere;
    }
    return Observation_t.NoMut;
  }
  
  private Observation_t emissionTypeHetInPair(string alleles, size_t ind1, size_t ind2) const {
    char dummy = '-';
    foreach(k, a; alleles) {
      if(k != ind1 && k != ind2) {
        if(dummy == '-')
          dummy = a;
        else {
          if(a != dummy)
            return Observation_t.DoubleMut;
        }
      }
    }
    return Observation_t.MutInPair;
  }
  
  double emissionProb(Observation_t flag, double t, double tTot) const {
    if(nrHaplotypes == 2) {
      if(flag == Observation_t.NoMut)
        return exp(-2.0 * mu * t);
      else
        return 1.0 - exp(-2.0 * mu * t);
    }
    
    if(to!double(nrHaplotypes) * t > tTot)
      // return 0.0;
      tTot = to!double(nrHaplotypes) * t;
    
    double first_coal_frac = 2.0 * t / tTot;
      
    double ret;
    final switch(flag) {
      case Observation_t.DoubleMut:
      ret = 0.0; // forbidden case
      break;
      case Observation_t.NoMut:
      ret = exp(-mu * tTot);
      break;
      case Observation_t.MutElsewhere:
      ret = (1.0 - exp(-mu * tTot)) * (1.0 - first_coal_frac) /
            (nrHaplotypes > 2 ? (2.0 ^^ (nrHaplotypes - 2) - 1.0) : 1.0);
      break;
      case Observation_t.MutInPair:
      ret = (1.0 - exp(-mu * tTot)) * first_coal_frac / 2.0;
      break;
    }
    return ret;
  }
}

unittest {
  writeln("test emissionRate.emissionType");
  auto emissionRate = new EmissionRate(4, 0.01);
  assert(emissionRate.emissionType("0000", 0, 1) == emissionRate.Observation_t.NoMut);
  assert(emissionRate.emissionType("0001", 0, 1) == emissionRate.Observation_t.MutElsewhere);
  assert(emissionRate.emissionType("0011", 0, 1) == emissionRate.Observation_t.MutElsewhere);
  assert(emissionRate.emissionType("0100", 0, 1) == emissionRate.Observation_t.MutInPair);
  assert(emissionRate.emissionType("1000", 0, 1) == emissionRate.Observation_t.MutInPair);
  assert(emissionRate.emissionType("1001", 0, 1) == emissionRate.Observation_t.DoubleMut);
}

unittest {
  writeln("test emissionRate.emissionProb");
  auto emissionRate = new EmissionRate(4, 0.01);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.DoubleMut, 1.0, 2.0) == 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.DoubleMut, 0.0, 2.0) == 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.DoubleMut, 0.5, 2.0) == 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.DoubleMut, 1.5, 2.0) == 0.0);

  assert(emissionRate.emissionProb(emissionRate.Observation_t.NoMut, 1.0, 2.0) > 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.NoMut, 1.0, 2.0) < 1.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.NoMut, 0.0, 2.0) > 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.NoMut, 0.0, 2.0) < 1.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.NoMut, 0.5, 2.0) > 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.NoMut, 0.5, 2.0) < 1.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.NoMut, 1.5, 2.0) > 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.NoMut, 1.5, 2.0) < 1.0);

  // assert(emissionRate.emissionProb(emissionRate.Observation_t.MutElsewhere, 1.0, 2.0) == 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.MutElsewhere, 0.0, 2.0) > 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.MutElsewhere, 0.0, 2.0) < 1.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.MutElsewhere, 0.5, 2.0) > 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.MutElsewhere, 0.5, 2.0) < 1.0);
  // assert(emissionRate.emissionProb(emissionRate.Observation_t.MutElsewhere, 1.5, 2.0) == 0.0);

  // assert(emissionRate.emissionProb(emissionRate.Observation_t.MutInPair, 1.0, 2.0) == 0.0);
  // assert(emissionRate.emissionProb(emissionRate.Observation_t.MutInPair, 0.0, 2.0) == 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.MutInPair, 0.5, 2.0) > 0.0);
  assert(emissionRate.emissionProb(emissionRate.Observation_t.MutInPair, 0.5, 2.0) < 1.0);
  // assert(emissionRate.emissionProb(emissionRate.Observation_t.MutInPair, 1.5, 2.0) == 0.0);
  
}
