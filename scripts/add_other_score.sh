#!/bin/bash

for file in data-sets/stimulated_vs_US28_*/proteins.nodes
do
./scripts/add_score.py -n ${file} -s deseq2/stimulated_vs_UL33_deseq2.txt 
done


for file in data-sets/stimulated_vs_UL33_*/proteins.nodes
do
./scripts/add_score.py -n ${file} -s deseq2/stimulated_vs_US28_deseq2.txt 
done
