#!/usr/bin/env python3

import argparse
import os
import sys
import random

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--nr_bootstraps", type=int, help="nr of bootstraps [20]", default=20)
parser.add_argument("-s", "--chunk_size", type=int, help="size of bootstrap chunks [5000000]", default=5000000)
parser.add_argument("--chunks_per_chromosome", type=int,
                    help="nr of chunks to put on one chromosome in the bootstrap [20]", default=20)
parser.add_argument("--nr_chromosomes", type=int, help="nr of chromosomes to write [30]", default=30)
parser.add_argument("--seed", type=int, help="initialize the random number generator")
parser.add_argument("out_dir_prefix", help="directory-prefix to write bootstraps to")
parser.add_argument("files", nargs="+")
args = parser.parse_args()

if args.seed is not None:
    random.seed(args.seed)

chunks = []
for fn in args.files:
    print("reading", fn, file=sys.stderr)
    chunks_in_chrom = []
    f = open(fn, "r")
    for line in f:
        fields = line.strip().split()
        pos, nr_called_sites = map(int, fields[1:3])
        alleles = fields[3]
        chunk_index = pos // args.chunk_size
        rel_pos = pos % args.chunk_size
        if nr_called_sites > rel_pos:
            nr_called_sites = rel_pos
        while chunk_index >= len(chunks_in_chrom):
            chunks_in_chrom.append([])
        chunks_in_chrom[chunk_index].append((rel_pos, nr_called_sites, alleles))
    chunks.extend(chunks_in_chrom)

print("created {} chunks".format(len(chunks)), file=sys.stderr)

for bootstrap_id in range(1, args.nr_bootstraps +1):
    for chr_ in range(1, args.nr_chromosomes + 1):
        chr_dir = "{}_{}".format(args.out_dir_prefix, bootstrap_id)
        if not os.path.exists(chr_dir):
            os.makedirs(chr_dir)
        chr_file = "{}/bootstrap_multihetsep.chr{}.txt".format(chr_dir, chr_)
        print("writing", chr_file, file=sys.stderr)
            
        f = open(chr_file, "wt")
        for i in range(args.chunks_per_chromosome):
            left_pos = i * args.chunk_size
            chunk_id = random.randrange(len(chunks))
            for pos, nr_called_sites, alleles in chunks[chunk_id]:
                print(chr_, left_pos + pos, nr_called_sites, alleles, sep="\t", file=f)
        f.close()
