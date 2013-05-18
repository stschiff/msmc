msmc_src = model/emission_rate.d model/transition_rate.d \
model/gsl_matrix_vector.d model/data.d model/time_intervals.d model/msmc_hmm.d model/propagation_core_naiveImpl.d \
model/propagation_core.d model/rate_integrator.d model/coalescence_rate.d \
model/triple_index_marginal.d model/triple_index.d model/msmc_model.d powell.d brent.d \
model/propagation_core_fastImpl.d maximization_step.d expectation_step.d baumwelch.d inference.d msmc.d \
../Utils/utils.d msmc_utils.d amoeba.d expectation.d maximization.d stats.d \
model/stateVec.d model/stateVecAllocator.d branchlength.d decode.d print.d

unittest_src = unittest.d model/emission_rate.d model/transition_rate.d \
model/gsl_matrix_vector.d model/data.d model/time_intervals.d model/msmc_hmm.d model/propagation_core_naiveImpl.d \
model/propagation_core.d model/rate_integrator.d model/coalescence_rate.d \
model/triple_index_marginal.d msmc_utils.d model/triple_index.d model/msmc_model.d powell.d brent.d \
model/propagation_core_fastImpl.d maximization_step.d expectation_step.d ../Utils/utils.d amoeba.d \
model/stateVec.d model/stateVecAllocator.d

memoryTest_src = model/emission_rate.d model/transition_rate.d \
model/gsl_matrix_vector.d model/data.d model/time_intervals.d model/msmc_hmm.d model/propagation_core_naiveImpl.d \
model/propagation_core.d model/rate_integrator.d model/coalescence_rate.d \
model/triple_index_marginal.d msmc_utils.d model/triple_index.d model/msmc_model.d powell.d brent.d \
model/propagation_core_fastImpl.d maximization_step.d expectation_step.d ../Utils/utils.d amoeba.d \
model/stateVec.d memoryTest.d model/stateVecAllocator.d

all : unittest build/msmc

unittest :
	dmd -unittest ${unittest_src} -L-lgsl -L-lgslcblas -L-L/opt/local/lib -odbuild -ofbuild/unittest
	build/unittest

build/msmc : ${msmc_src}
	dmd ${msmc_src} -O -L-lgsl -L-lgslcblas -L-L/opt/local/lib -odbuild -ofbuild/msmc

build/memoryTest : ${memoryTest_src}
	dmd ${memoryTest_src} -L-lgsl -L-lgslcblas -L-L/opt/local/lib -odbuild -ofbuild/memoryTest

clean :
	rm build/*
