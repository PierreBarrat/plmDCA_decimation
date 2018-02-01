#!/bin/bash

steps=(0 5 10 15 20 25 30 35 40)

while read f; do
	echo ""
	echo "Family: $f"
	cd $f"_results"
	mkdir -p samples
	lseq=`head -n 1 msa_numerical.txt |awk '{print NF}'`
	echo "Length of sequence: $lseq"
	cd samples
	for s in ${steps[*]}; do
		echo "Decimation: $s"
		st="$(date -u +%s)" 
		../../sources/graph_mat_to_MCMC ../decimation_results/plmInf_$s"_mat.txt" plmInf_$s"_MCMC.txt" $lseq 21
		../../sources/do_montecarlo -n $lseq -q 21 -t 5000 -m 5000 -s $s < plmInf_$s"_MCMC.txt" &> /dev/null
		rm out_energies* 
		mv out_samp* sample_$s.txt
		et="$(date -u +%s)"
		elapsed="$(($et-$st))"
		rm plmInf_$s"_MCMC.txt" 
		echo "$elapsed seconds"
	done
	echo "${steps[*]}" > steps.txt
	echo ""
done <$1