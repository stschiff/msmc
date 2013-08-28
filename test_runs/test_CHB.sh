#!/usr/bin/env bash

THETA=0.00069863
MU=0.0003493
INPUTFILES=$(ls $lhome/CG69/hetsep_combined/CHB_combined/CHB_combined_multihetsep_phased_chr*.txt | tr '\n' ' ')
echo "msmc -m $MU -o tmp/CHB_test_run -t 8 -R $INPUTFILES"