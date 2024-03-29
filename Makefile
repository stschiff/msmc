# Set this variable to your static gsl libraries
# GSL = /usr/lib/libgsl.a /usr/lib/libgslcblas.a
GSL=/usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a

all: build/msmc build/decode

build/msmc : model/*.d powell.d brent.d maximization_step.d expectation_step.d msmc.d branchlength.d logger.d
	dmd -debug -O ${GSL} -odbuild -ofbuild/msmc $^

build/maximize : model/*.d powell.d brent.d maximization_step.d logger.d maximize.d
	dmd -O ${GSL} -odbuild -ofbuild/maximize $^

build/test/msmc : model/*.d powell.d brent.d maximization_step.d expectation_step.d msmc.d branchlength.d logger.d
	dmd -O ${GSL} -odbuild/test -ofbuild/test/msmc $^

build/decode : model/*.d decode.d branchlength.d
	dmd ${GSL} -odbuild -ofbuild/decode $^

.PHONY : unittest testcoverage clean

testcoverage : model/*.d unittest.d powell.d brent.d maximization_step.d expectation_step.d amoeba.d logger.d
	mkdir -p code_coverage
	dmd -unittest -cov ${GSL} -odbuild -ofbuild/unittest $^
	build/unittest
	mv *.lst code_coverage/

unittest : model/*.d unittest.d powell.d brent.d maximization_step.d expectation_step.d logger.d branchlength.d
	dmd -debug -unittest ${GSL} -odbuild -ofbuild/unittest $^
	build/unittest

clean :
	find build -type f -delete

