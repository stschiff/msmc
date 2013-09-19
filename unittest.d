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
import model.emission_rate;
import model.transition_rate;
import model.gsl_matrix_vector;
import model.data;
import model.time_intervals;
import model.msmc_hmm;
import model.propagation_core_naiveImpl;
import model.propagation_core;
import model.rate_integrator;
import model.coalescence_rate;
import model.triple_index_marginal;
import model.triple_index;
import model.msmc_model;
import model.propagation_core_fastImpl;
import model.stateVec;
import model.stateVecAllocator;
import powell;
import brent;
import maximization_step;
import expectation_step;
import amoeba;
import branchlength;

void main() {
  writeln("all tests processed");
}