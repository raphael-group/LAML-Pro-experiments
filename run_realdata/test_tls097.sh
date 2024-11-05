#!/bin/bash
#

input_observations="/Users/gc3045/laml2_experiments/proc_realdata/TLS_097_test.json"
topology="random_topology.nwk"
priors="/Users/gc3045/laml2_experiments/tlsdata/TLS_097_unfiltered_cm_priors.pickle"
priors="/Users/gc3045/laml2_experiments/proc_realdata/TLS_097_priors.json"
#priors="/Users/gc3045/laml2_experiments/tlsdata/TLS_097_unfiltered_cm_priors.pickle"
outputfile="LAML2_test"

echo "python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p $priors -m -1 -y "allele_counts" -M "PMMC" -o ${outputfile} -v"
python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p $priors -m -1 -y "allele_counts" -M "PMMC" -o ${outputfile} -v
