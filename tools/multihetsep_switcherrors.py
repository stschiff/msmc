#!/usr/bin/env python3

import argparse
import random

def flip_phase(alleles, ind):
    new_alleles = ""
    nr_ind = len(alleles) // 2
    for i in range(nr_ind):
        if i == ind:
            new_alleles += alleles[i * 2 + 1] + alleles[i * 2]
        else:
            new_alleles += alleles[i * 2] + alleles[i * 2 + 1]
    return new_alleles

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seed", type=int, help="Random Seed")
parser.add_argument("length", type=int, help="Expected length (in bp) of correctly phased segments")
parser.add_argument("file", help="Input file")
args = parser.parse_args()

if args.seed is not None:
    random.seed(args.seed)

f = open(args.file, "r")

rate = 1.0 / args.length

nr_ind = None
first_pos = None
next_switch = []
switch_state = []
for line in f:
    fields = line.strip().split()
    chrom = fields[0]
    pos = int(fields[1])
    nr_called_sites = int(fields[2])
    alleles = fields[3]

    if nr_ind is None:
        nr_ind = len(alleles) // 2
        first_pos = pos
        for i in range(nr_ind):
            next_switch.append(int(random.expovariate(rate)))
            switch_state.append(False)
    for i in range(nr_ind):
        if pos - first_pos >= next_switch[i]:
            switch_state[i] = not switch_state[i]
            next_switch[i] += int(random.expovariate(rate))
        if switch_state[i]:
            alleles = flip_phase(alleles, i)
    
    # print(chrom, pos, nr_called_sites, alleles, fields[3], switch_state, sep="\t")
    print(chrom, pos, nr_called_sites, alleles, sep="\t")