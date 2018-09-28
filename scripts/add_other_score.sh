#!/bin/bash

for file in data-sets/US28_vs_stimulated_*/proteins.nodes
do
./scripts/add_score.py -n ${file} -s deseq2/UL33_vs_stimulated_deseq2.txt 
done


for file in data-sets/UL33_vs_stimulated_*/proteins.nodes
do
./scripts/add_score.py -n ${file} -s deseq2/US28_vs_stimulated_deseq2.txt 
done
