# The multiple sequentially Markovian coalescent (MSMC)

This software implements MSMC, a method to infer population size and gene flow from multiple genome sequences. The accompanying manuscript (Schiffels and Durbin, 2013) that describes the method and its analysis to human genome sequences is currently under review.

In short, msmc can infer

* the scaled population size of a single population
* gene flow across multiple populations

as a function of time from multiple phased haplotypes. When only two haplotypes are given, MSMC is similar to [PSMC](http://github.com/lh3/psmc), and we call it PSMC' because of subtle differences in the method and the underlying model, which allows PSMC' to infer more accurately the recombination rate.

MSMC infers the lowest part of the genealogical tree, up to the first coalescence, at every position along the genome. Therefore, the number of samples strongly affects the time of inference: if four phased haplotypes are analyzed (e.g. from two diploid samples), the typical time of first coalescence is 6 times more recent than with PSMC. With 8 haplotypes, the time is already 28 times more recent. Therefore, in contrast to many other methods, the rule "more samples => higher resolution" is not true here. Instead, it is recommended to run MSMC on multiple number of samples, beginning with 2, then increasing the sample size up to 8, if available. 

We tested MSMC on up to 8 haplotypes. The computational complexity of the model so far limits its application to possibly a dozen of haplotypes.

# Compilation

To build and run MSMC, the [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/) must be installed on your system. (This dependence is likely to be dropped at some point soon).

To build the program on a mac, simply run `make -f Makefile.mac` in the directory in which your copy of the source code is located. The program is written in the [D programming language](http://dlang.org). The reference compiler from Digitalmars can be downloaded [here](http://dlang.org/download.html).

A precompiled version will be made available for Mac OS X and for Linux. Please check later.

# Usage for population size estimates
## Input files

MSMC takes as input several files, one for each chromosome, each with a list of segregating sites, including the number of homozygous called sites between. Here is an example bit of an input file for MSMC:

    1 58432 63 TCCC
    1 58448 16 GAAA
    1 68306 15 CTTT
    1 68316 10 TCCC
    1 69552 8 GCCC
    1 69569 17 TCCC
    1 801848 9730 CCCA
    1 809876 1430 AAAG
    1 825207 1971 TCCT,CTTC
    1 833223 923 TCCC

The four columns are:

1. the chromosome (can be any string)
2. the position in the chromosome
3. the number of called homozygous sites since the last segregating site, including the given location
4. the ordered and phased pattern of alleles at multiple haplotypes. Multiple observations can be given, separated by a comma to indicate ambiguous haplotype phasing. Unknown alleles are indicated by "?", e.g. for missing data, or for sites which cannot be phased.

Input files can be generated from BAM files, VCF files and other sources relatively straight forward. The necessary scripts will be made available soon.

For the following, we need an estimate of the diversity from the input files. A standard estimator for θ(=4Nµ), along with several other statistics can be obtained with the command

    msmc stats <file1> <file2> ...

## Estimating the total branchlength
This step is necessary if more than two haplotypes are analyzed. Local estimates for the scaled total branchlength of local genealogies along the sequences can be obtained with the following command

    msmc branchlength -m <mutation_rate> <input_file> > out_file

which must be run for every input file. The scaled mutation rate. 2Nµ should be set to half the estimate obtained in the previous step for θ. The output files look like this:

    1 58448 16 GAAA 0.0766465
    1 68306 15 CTTT 8.58025
    1 68316 10 TCCC 8.58025
    1 69552 8 GCCC 8.58025
    1 69569 17 TCCC 8.58025
    1 801848 9730 CCCA 0.573957
    1 809876 1430 AAAG 0.84195
    1 825207 1971 TCCT,CTTC 2.19634
    1 833223 923 TCCC 5.78576
    1 833824 100 CTTT 5.78576

where the local estimate of the total branchlength is added to every segregating site as fifth column.

## Main inference
We can now run the main program, `msmc inference`, on our annotated datafile (or in case of 2 haplotypes without the branchlength annotation). To run the program on a machine with 8 processor cores, the minimal command line looks like this:

    msmc inference -m <mutation_rate> -t 8 --fixedRecombination -o <out_file> <input_file1> <input_file2> ...

It is very important that the mutation rate is the _same_ as was used in the previous step for `msmc branchlength`!

## Printing the output
The output from `msmc inference` is in JSON format, which is a simple key-value format that is human readable. A helper program to print in a simple tabular format is available as `msmc print` and can be called like this:

    msmc print -s -y <out_file>

The first lines of the tabular output look like this:

    0	1.92217e+06
    3623.48	138619
    7341.09	22498.9
    11157.9	15691.5
    15079.2	13743.6
    19111	11315.2
    23259.7	9071.03
    27532.3	7351.55
    31936.3	6154.67
    36480.1	5354.62

where the first column contains the left boundary of the time intervals, and the second column the population size estimates.

#Usage for gene flow estimates
... to be continued ...
