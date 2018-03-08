#Load required modules
import argparse
import os
import pandas as pd
import scipy.stats as stats

parser = argparse.ArgumentParser(description = 'Computes KEGG pathway enrichment')
parser.add_argument('-m', type = str, help = 'module file', required=True)
parser.add_argument('-b', type = str, help = 'background file', required=True)


module = set()
background = set()

try:
    args = parser.parse_args()
except IOError:
    parser.error()

with open(args.m, 'rt') as nodes:
    for line in nodes:
        module.add(line.split()[0])

with open(args.b, 'rt') as nodes:
    for line in nodes:
        if not line.startswith('#'): background.add(line.split()[0])

df = pd.read_csv('KEGG_pathways.csv')

for index, row in df.iterrows():
    pathway = set(eval(row['genes']))
    pathway = pathway & background
    module_pathway = len(module & pathway)
    module_non_pathway = len(module - pathway)
    non_module_pathway = len(pathway - module)
    non_module_non_pathway = len(background - (module | pathway))

    oddsratio, pvalue = stats.fisher_exact([[module_pathway, non_module_pathway], [module_non_pathway, non_module_non_pathway]])
    print(pvalue, module_pathway, module_non_pathway, non_module_pathway, non_module_non_pathway, row['description'])