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

module model.triple_index;
unittest {
  auto sI = new TripleIndex(4, 4);
  assert(sI.getIndexFromTriple(Triple(0, 0, 1)) == 0);
  assert(sI.getIndexFromTriple(Triple(1, 0, 2)) == 7);
  assert(sI.getIndexFromTriple(Triple(3, 2, 3)) == 23);
}

unittest {
  auto sI = new TripleIndex(4, 4);
  assert(sI.getTripleFromIndex(0) == Triple(0, 0, 1));
  assert(sI.getTripleFromIndex(7) == Triple(1, 0, 2));
  assert(sI.getTripleFromIndex(23) == Triple(3, 2, 3));
}

struct Triple {
  size_t time;
  size_t ind1;
  size_t ind2;
}

class TripleIndex {
  size_t nrTimeIntervals, nrIndividuals;
  
  private size_t[][][] tripleToIndexMap;
  private Triple[] indexToTripleMap;
  
  this(size_t nrTimeIntervals, size_t nrIndividuals) {
    this.nrTimeIntervals = nrTimeIntervals;
    this.nrIndividuals = nrIndividuals;
    
    buildMaps();
  }
  
  @property size_t nrStates() const {
    return nrTimeIntervals * nrIndividuals * (nrIndividuals - 1) / 2;
  }
  
  private void buildMaps() {
    tripleToIndexMap = new size_t[][][](nrTimeIntervals, nrIndividuals, nrIndividuals);
    indexToTripleMap = new Triple[nrStates()];
    auto index = 0;
    foreach(t; 0 .. nrTimeIntervals) foreach(i; 0 .. nrIndividuals - 1) foreach(j; i + 1 .. nrIndividuals) {
      tripleToIndexMap[t][i][j] = index;
      indexToTripleMap[index] = Triple(t, i, j);
      ++index;
    }
  }
  
  size_t getIndexFromTriple(Triple triple) const
    in {
      assert(triple.ind1 < triple.ind2);
      assert(triple.time < nrTimeIntervals);
    }
    out(index) {
      assert(index < nrStates());
    }
  do {
    return tripleToIndexMap[triple.time][triple.ind1][triple.ind2];
  }
  
  Triple getTripleFromIndex(size_t index) const
    in {
      assert(index < nrStates());
    }
    out(triple) {
      assert(triple.ind1 < triple.ind2);
      assert(triple.time < nrTimeIntervals);
    }
  do {
    return indexToTripleMap[index];
  }
}
