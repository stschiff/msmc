#!/usr/bin/env bash

THETA=0.00069863
MU=0.0003493
# DIR=$lhome/CG69/hetsep_combined
DIR=$HOME/Data
INPUTFILES=($(ls $DIR/CHB_combined/CHB_combined_multihetsep_phased_chr*.txt | tr '\n' ' '))
echo "msmc -m $MU -o tmp/CHB_test_run -t 2 -R -p 20*1 ${INPUTFILES[@]:0:8}"