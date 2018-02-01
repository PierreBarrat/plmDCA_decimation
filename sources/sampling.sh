#!/bin/bash

steps=(0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95)
for s in ${steps[*]}; do
	echo $s
done

while read f; do
	echo $f
	cd $f_results
	mkdir -p samples
	lseq=`head -n 1 msa_numerical.txt |awk '{print NF}'`
	echo "Length of sequence: $lseq"
	cd samples
	for s in ${steps[*]}; do
		../../sources/graph_mat_to_MCMC ../decimation_results/plmInf_$s_mat.txt plmInf_$s_MCMC.txt $lseq 21
		../../sources/do_montecarlo -n $lseq -q 21 -t 5000 -m 2000 < plmInf_$s_MCMC.txt 
		rm out_energies* 
		mv out_samp* sample_$s.txt
	done
done <$1