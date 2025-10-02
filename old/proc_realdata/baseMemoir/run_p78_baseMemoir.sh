#!/bin/bash

topology="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/p78_starting_tree.nwk"
input_observations="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/baseMemoir.msa.txt"
emission_file="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/baseMemoir.emissions.txt"
outputfile="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/baseMemoir_runs/LAML2_p78_baseMem"

python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p "uniform" -m -1 -y "character_matrix" -M "PMMN" -o ${outputfile} -v --noSilence --e $emission_file --topology_search

