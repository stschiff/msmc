msmc_src = model/emission_rate.d model/transition_rate.d \
model/gsl_matrix_vector.d model/data.d model/time_intervals.d model/msmc_hmm.d model/propagation_core_naiveImpl.d \
model/propagation_core.d model/rate_integrator.d model/coalescence_rate.d \
model/triple_index_marginal.d model/triple_index.d model/msmc_model.d powell.d brent.d \
model/propagation_core_fastImpl.d maximization_step.d expectation_step.d baumwelch.d inference.d msmc.d \
utils.d amoeba.d expectation.d maximization.d stats.d \
model/stateVec.d model/stateVecAllocator.d branchlength.d decode.d print.d

unittest_src = unittest.d model/emission_rate.d model/transition_rate.d \
model/gsl_matrix_vector.d model/data.d model/time_intervals.d model/msmc_hmm.d model/propagation_core_naiveImpl.d \
model/propagation_core.d model/rate_integrator.d model/coalescence_rate.d \
model/triple_index_marginal.d model/triple_index.d model/msmc_model.d powell.d brent.d \
model/propagation_core_fastImpl.d maximization_step.d expectation_step.d utils.d amoeba.d \
model/stateVec.d model/stateVecAllocator.d

all : unittest debug

unittest :
	dmd -unittest ${unittest_src} -L-lgsl -L-lgslcblas -L-L/opt/local/lib -odbuild -ofbuild/unittest
	build/unittest

debug : build/debug/msmc

release : build/release/msmc

build/debug/msmc : ${msmc_src}
	dmd ${msmc_src} -L-lgsl -L-lgslcblas -L-L/opt/local/lib -odbuild -ofbuild/debug/msmc

build/release/msmc : ${msmc_src}
	ldmd2 ${msmc_src} -O -release -inline -L-lgsl -L-lgslcblas -L-L/opt/local/lib -odbuild -ofbuild/release/msmc

testcoverage :
	mkdir -p code_coverage
	dmd -unittest -cov ${unittest_src} -L-lgsl -L-lgslcblas -L-L/opt/local/lib -odbuild -ofbuild/unittest
	build/unittest
	mv *.lst code_coverage/

clean :
	find build -type f -delete

.PHONY : all unittest testcoverage clean