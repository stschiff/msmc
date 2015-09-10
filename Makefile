# Set this variable to your static gsl libraries
GSL = /usr/lib/libgsl.a /usr/lib/libgslcblas.a
# GSL = /opt/local/lib/libgsl.a /opt/local/lib/libgslcblas.a

build/msmc : model/*.d powell.d brent.d maximization_step.d expectation_step.d msmc.d branchlength.d logger.d
	dmd -O ${GSL} -odbuild -ofbuild/msmc $^

build/maximize : model/*.d powell.d brent.d maximization_step.d logger.d maximize.d
	dmd -O ${GSL} -odbuild -ofbuild/maximize $^

build/test/msmc : model/*.d powell.d brent.d maximization_step.d expectation_step.d msmc.d branchlength.d logger.d
	dmd -O ${GSL} -odbuild/test -ofbuild/test/msmc $^

build/test : model/*.d test.d
	dmd ${GSL} -odbuild -ofbuild/test $^

build/decode : model/*.d decode.d branchlength.d
	dmd ${GSL} -odbuild -ofbuild/decode $^

build/decodeStates : model/*.d decodeStates.d branchlength.d
	dmd ${GSL} -odbuild -ofbuild/decodeStates $^

.PHONY : unittest testcoverage clean

testcoverage : model/*.d unittest.d powell.d brent.d maximization_step.d expectation_step.d amoeba.d logger.d
	mkdir -p code_coverage
	dmd -unittest -cov ${GSL} -odbuild -ofbuild/unittest $^
	build/unittest
	mv *.lst code_coverage/

unittest : model/*.d unittest.d powell.d brent.d maximization_step.d expectation_step.d logger.d branchlength.d
	dmd -unittest ${GSL} -odbuild -ofbuild/unittest $^
	build/unittest

clean :
	find build -type f -delete

