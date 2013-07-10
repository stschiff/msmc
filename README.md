# The multiple sequentially Markovian coalescent (MSMC)

This software implements MSMC, a method to infer population size and gene flow from multiple genome sequences (Schiffels and Durbin, submitted).

In short, msmc can infer

* the scaled population size of a single population
* gene flow across multiple populations

as a function of time from multiple phased haplotypes. When only two haplotypes are given, MSMC is similar to [PSMC](http://github.com/lh3/psmc), and we call it PSMC' because of subtle differences in the method and the underlying model, which allows PSMC' to infer more accurately the recombination rate.

# Installation

Precompiled versions for Mac and Linux (both 64 bit) can be downloaded via ftp from

    ftp://ftp.sanger.ac.uk/pub/users/ss27/msmc/

To build MSMC yourself, the [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/) must be installed on your system.

To build the program, have a look at the two Makefiles. Adjust the path to the GSL and eventually run the `release` target. The program is written in the [D programming language](http://dlang.org). The reference compiler from Digitalmars can be downloaded [here](http://dlang.org/download.html).

The main program can be called via the command `msmc`. It outputs several subcommands when called without further options. Each subprogram outputs all its options if called without options.

# Preparing the input files

MSMC takes as input several files, one for each chromosome, each with a list of segregating sites, including the number of homozygous called sites between. Here is an example bit of an input file for MSMC:

    1   58432	63	TCCC
    1   58448	16	GAAA
    1	68306	15	CTTT
    1	68316	10	TCCC
    1	69552	8	GCCC
    1	69569	17	TCCC
    1	801848	9730	CCCA
    1	809876	1430	AAAG
    1	825207	1971	TCCT,CTTC
    1	833223	923	TCCC

The four (tab-separated) columns are:

1. the chromosome (can be any string)
2. the position in the chromosome
3. the number of called homozygous sites since the last segregating site, including the given location
4. the ordered and phased pattern of alleles at multiple haplotypes. Multiple observations can be given, separated by a comma to indicate ambiguous haplotype phasing. Unknown alleles are indicated by "?", e.g. for missing data.

Input files can be generated from BAM files, VCF files and other sources relatively straight forward. Useful scripts can be found in the `tools` directory (in order to run them, the D-compiler must be installed and rdmd somewhere in the path). 

# Estimating the total branchlength

For the following we need an estimate of the diversity from the input files. A standard estimator for θ(=4Nµ), along with several other statistics can be obtained with the command

    msmc stats <file1> <file2> ...

If more than two haplotypes are analyzed, local estimates for the scaled total branchlength of local genealogies along the sequences must be obtained with the following command

    msmc branchlength -m <mutation_rate> <input_file> > out_file

which must be run for every input file. The scaled mutation rate, 2Nµ, should be set to half the estimate obtained in the previous step for θ. The output files look like this:

    1   58448	16	GAAA	0.0766465
    1	68306	15	CTTT	8.58025
    1	68316	10	TCCC	8.58025
    1	69552	8	GCCC	8.58025
    1	69569	17	TCCC	8.58025
    1	801848	9730	CCCA	0.573957
    1	809876	1430	AAAG	0.84195
    1	825207	1971	TCCT,CTTC	2.19634
    1	833223	923	TCCC	5.78576
    1	833824	100	CTTT	5.78576

where the local estimate of the total branchlength is added to every segregating site as fifth column.

# Estimation of Historical Population sizes

We can now run the main program, `msmc inference`, on our annotated datafile (or in case of 2 haplotypes without the branchlength annotation). To run the program on a machine with 8 processor cores, the minimal command line looks like this:

    msmc inference -m <mutation_rate> -t 8 --fixedRecombination -o <out_file> <input_file1> <input_file2> ...

It is very important that the mutation rate is the _same_ as was used in the previous step for `msmc branchlength`!

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

# Estimation of historical gene flow

In this case, the data file should consist of the alleles from two subpopulations, e.g. 4 haplotypes with two haplotypes from each subpopulation. The command for running the inference is then

    msmc inference -m <mutation_rate> -t 8 --fixedRecombination -P 0,0,1,1 -o <out_file> <input_file1> <input_file2> ...

where the flag `-P 0,0,1,1` specifies that the four alleles are sampled from two subpopulations `0` and `1`.

To output the gene flow estimates directly, use

    msmc print -s -y -w crossLambda <out_file>

