# Short guide to MSMC

## Input File Format

MSMC takes as input several files, one for each chromosome, each with a list of segregating sites, including a column to denote how many sites have been called since the last segregating site. Here is an example bit of an input file for MSMC:

    1	58432	63	TCCC
    1	58448	16	GAAA
    1	68306	15	CTTT
    1	68316	10	TCCC
    1	69552	8	GCCC
    1	69569	17	TCCC
    1	801848	9730	CCCA
    1	809876	1430	AAAG
    1	825207	1971	CCCT,CCTC
    1	833223	923	TCCC

The four (tab-separated) columns are:

1. the chromosome (can be any string)
2. the position in the chromosome
3. the number of called sites (homozygous, except the site itself which can be hom. or het.) since the last segregating site. This number *includes* the given location. This means that this number must always be greater than zero! It also means that this number cannot be larger than the difference between the current and the previous position.
4. the ordered and phased alleles of the multiple haplotypes. If phasing is unknown, multiple phasings can be given, separated by a comma to indicate the different possibilities. This is the case in the line before the last in the example above. Unknown alleles can be indicated by "?", but they can also simply be left out and expressed through a reduced number of called sites in the line of the next variant.

## Generate input files

Generating input files for msmc normally consists of these steps:

1. Starting from either a bam-file or a complete Genomics masterVarBeta file, you generate a single-sample VCF and a mask-file. This can be done by using the scripts `bamCaller.py` or `cgCaller.py` in the [msmc-tools](http://github.com/stschiff/msmc-tools)-Repository.
2. If your samples are unrelated and you want to run MSMC on more than 2 haplotypes at a time, you need to phase the VCFs using a tools like `shapeit`. You find a script called `run_shapeit.sh` in [msmc-tools](http://github.com/stschiff/msmc-tools), but you will need to customize it and consult the shapeit documentation to understand it and run it successfully. When you have trios, you can omit this step. Note that you can also decide to run MSMC on more than 2 unphased haplotypes, which has been shown to generate some biases, in particular when you study population divergences via the relative gene flow analysis as described below. Please check out Schiffels and Durbin, 2014.
3. Generate the input file using `generate_multihetsep.py` in [msmc-tools](http://github.com/stschiff/msmc-tools).

## Estimation of historical effective population sizes

To run the program, the minimal command line looks like this:

    msmc --fixedRecombination -o my_msmc_output file1.txt file2.txt file3.txt [...]

For running on only one individual (two haplotypes), you should leave out the flag `--fixedRecombination`. If no mutation rate is given, as in the line above, Watterson's estimator is used to determine theta. If you find that confusing, don't worry about it, it simply sets the exact placement of the time-intervals. The flag "--fixedRecombination" is recommended, but population size estimates are pretty robust whether or not this flag is set.
More command line options are printed when simply typing `msmc`.

The program outputs three files:
* A file called something.log: This file contains the same logging information that is printed out while it runs.
* A file called something.loop.txt. This file contains a table with parameter estimates after each iteration step. The columns of the table are: the recombination rate (fixed in the above example), the log-likelihood, and a comma-separated list of coalescence rates.
* A file called something.final.txt. This file contains a table with the final parameter estimates.

The final file contains multiple columns with a header line. Here are the first rows of an example:


    time_index	left_time_boundary	right_time_boundary	lambda_00
    0	-0	2.09028e-06	1086.3
    1	2.09028e-06	4.23486e-06	3373.81
    2	4.23486e-06	6.43663e-06	3726.96
    3	6.43663e-06	8.69874e-06	3009.26

The first column simply gives the zero-based index of the time interval. The second and third columns give the scaled begin and end time of the interval. The fourth column gives the scaled inverse population size (scaled coalescence rate) of the interval. See below on how to convert scaled times and rates to real numbers.

## Analyzing population separations

If the dataset contains samples from different subpopulations of the same species (i.e. not too diverged), MSMC can estimate gene flow as a function of time. Assuming you have four haplotypes, two of which are from population 0, and two from population 1, the command for running the inference is

    msmc --fixedRecombination --skipAmbiguous -P 0,0,1,1 -o my_msmc_output file1.txt file2.txt file3.txt [...]

Here, the flag `-P 0,0,1,1` specifies that the four alleles are sampled from two subpopulations `0` and `1`. These need to be given as numbers starting from 0. For eight haplotypes, 4 sampled from each subpopulation, this would read `-P 0,0,0,0,1,1,1,1`. MSMC has been tested only on an equal number of haplotypes in each of the two subpopulation so far, feel free to try unequal numbers. In principle, MSMC allows more than two subpopulations as well, but that hasn't been tested yet.

Note that population separation analysis is expensive in memory and time. You might want to reduce the resolution for some test runs to 30 segments (using e.g. `-p 8*1+11*2`), or even 20 segments (using e.g. `-p 20*1`).

The flag `--skipAmbiguous` is recommended for gene flow estimation. It means that sites with ambiguous phasing, as described above, are removed from the analysis.

The output file now contains several coalescence rate estimates. Here is an example bit:

    time_index	left_time_boundary	right_time_boundary	lambda_00	lambda_01	lambda_11
    0	-0	2.79218e-06	2605.47	71.9887	4206.61
    1	2.79218e-06	5.68236e-06	6451.92	1256.07	3897.26
    2	5.68236e-06	8.67766e-06	3152.31	736.499	2790.45
    3	8.67766e-06	1.1786e-05	2526.36	1075.56	2790.33

Now the three columns titled lambda_?? denote the coalescence rates within and across the subpopulations. To get relative gene flow, you can compute the relative cross-coalescence rate: 2 * lambda01 / (lambda00 + lambda11).

## Scaling to real time and population sizes

MSMC outputs times and rates scaled by the mutation rate per basepair per generation.
First, scaled times are given in units of the per-generation mutation rate. This means that in order to convert scaled times to generations, divide them by the mutation rate. In humans, we used mu=1.25e-8 per basepair per generation.To convert generations into years, multiply by the generation time, for which we used 30 years.

To get population sizes out of coalescence rates, first take the inverse of the coalescence rate, scaledPopSize = 1 / lambda00. Then divide this scaled population size by 2*mu (yes, this factor 2 is different from the time scaling, sorry).

## Estimating the local tMRCA states

MSMC contains a tool called `decode`. It can be compiled by running `make build/decode`, which will put the executable under the `build` directory. This program takes as input file the standard MSMC input file for a single chromosome and outputs the posterior probabilities for each state. You run it like this:

    build/decode -m 0.0005 -r 0.0004 <input_file> > posteriors.txt

Type `build/decode` for a short help. The options `-m` and `-r` specify the scaled mutation- and recombination rates, respectively. Ideally, the scaled mutation rate should be set to half the heterozygosity in your data (half because here we scale with 2N, and the heterozygosity is scaled by 4N). In practice, for this tool the precise value is not too important. For example, in humans, I usually just assume a heterozygosity of 0.001 for Africans, and 0.00075 for Non-Africans, so I put either `-m 0.0005` for Africans or `-m 0.000375` for Non-Africans. The recombination rate does not affect results strongly, either, and I usually put it to to about 80% of the mutation rate, so `-r 0.0004` for Africans and `-r 0.0003` or so for Non-Africans.

The first row of the output is a header row specifying the lower ends of the time boundaries for each state. The rest of the output is a large matrix. At every row, the first value is the position in the genome, and the values from the second position are the posterior probabilities for all hidden tMRCA states in the model. The positions are by default spaced every 1000 basepairs, but that can be controlled with the `-s` flag.

Note that the number of states only corresponds to the number of time intervals if your input file contains only two haplotypes (PSMC mode). For more than 2 haplotypes, MSMC's hidden states also span the possible pairs of first coalescence. So for example, for four haplotypes, every time interval contains 6 hidden states, corresponding to the pairings 0,1; 0,2; 0,3; 1,2; 1,3; 2,3; where 0 denotes the first haplotype and 3 denotes the fourth haplotype. So with 4 time intervals (default) you have 240 hidden states, all of which are output in this matrix. Usually, I think this tool makes most sense in PSMC mode, so with just two input haplotypes.

As a first step, I recommend plotting the matrix as an image plot with squares for each matrix entry. This would then space all the states equally. It is more difficult to make a plot with the correct time boundaries, as that requires plotting squares with variable height (this is what has been done in Figure 1b of the paper).
