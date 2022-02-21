# The multiple sequentially Markovian coalescent (MSMC)

This software implements MSMC, a method to infer population size and gene flow from multiple genome sequences ([Schiffels and Durbin, 2014, Nature Genetics](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3015.html), or [Preprint](http://biorxiv.org/content/early/2014/05/21/005348)).

In short, msmc can infer

* the scaled population size of a single population as a function of time
* the timing and nature of population separations between two populations

from multiple phased haplotypes. When only two haplotypes are given, MSMC is similar to [PSMC](http://github.com/lh3/psmc), and we call it PSMC' because of subtle differences in the method and the underlying model, which allows PSMC' to infer more accurately the recombination rate.

# Installation and Requirements

Precompiled versions for Mac and Linux (both 64 bit) can be downloaded on the "Releases" tab within this github repository.

To build MSMC yourself, the [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/) must be installed on your system.

To build the program, have a look at the two Makefiles. Adjust the path to the GSL and eventually run `make -f Makefile.linux` or `make -f Makefile.mac`, respectively. The program is written in the [D programming language](http://dlang.org). The reference compiler from Digitalmars can be downloaded [here](http://dlang.org/download.html).

For generating the input files using my scripts, you need Python 3.4. I am sorry for this cutting edge dependency, I may make things compatible with Python 3.2 soon, but at the moment apparently my scripts won't work unless you use python 3.4.

# Getting Help

A general guide can be found [here](https://github.com/stschiff/msmc/blob/master/guide.md)

To get help, please submit an [issue](https://github.com/stschiff/msmc/issues), or [email me directly](https://www.eva.mpg.de/archaeogenetics/staff/stephan-schiffels/)
