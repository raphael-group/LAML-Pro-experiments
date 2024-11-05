#!/bin/bash

topology="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/baseMemoir_runs/LAML2_p910_baseMem_trees.nwk"
input_observations="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/baseMemoir.msa.txt"
emission_file="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/baseMemoir.emissions.txt"
outputfile="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/impute_sequences/LAML2_p910_baseMem_computellh"

python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p "uniform" -m -1 -y "character_matrix" -M "PMMN" -o ${outputfile} -v --noSilence --e $emission_file --topology_search --compute_llh --fixedParams "lambda=0.5011200336299421 phi=0.33964646464646403 nu=0 rho=0.9627415744955495"

