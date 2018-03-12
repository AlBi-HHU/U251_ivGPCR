#Load required modules
import argparse
import os
import pandas as pd
import scipy.stats as stats

parser = argparse.ArgumentParser(description = 'Computes KEGG pathway enrichment')
parser.add_argument('-m', type = str, help = 'module file', required=True)
parser.add_argument('-b', type = str, help = 'background file', required=True)
parser.add_argument('-ao', type = str, help = 'annotations output file', required=True)
parser.add_argument('-mo', type = str, help = 'memberships output file', required=True)

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

P = pd.DataFrame(columns = ['p', 'ID', 'description', 'genes']) # pathways enriched in module

for index, row in df.iterrows():
    pathway = set(eval(row['genes']))
    pathway = pathway & background
    module_pathway = len(module & pathway)
    module_non_pathway = len(module - pathway)
    non_module_pathway = len(pathway - module)
    non_module_non_pathway = len(background - (module | pathway))

    oddsratio, pvalue = stats.fisher_exact([[module_pathway, non_module_pathway], [module_non_pathway, non_module_non_pathway]])
    #P.append(pvalue, row['ID'], row['description'], list(module & pathway))
    P.at[index, 'p'] = pvalue
    P.at[index, 'ID'] = row['ID']
    P.at[index, 'description'] = row['description'][:-23]
    P.at[index, 'genes'] = list(module & pathway)

    #print(pvalue, module_pathway, module_non_pathway, non_module_pathway, non_module_non_pathway, row['description'])

P = P.sort_values(by=['p'])

with open (args.ao, 'at') as annotation:
    i = 0
    for index, row in P.iterrows():
        if i >= 30: break # only max. 30 pathways
        annotation.write(row['ID'] + "\tpathway\t" + str(row['p']) + "\t" + row['description'] + "\n")
        i = i + 1

with open (args.mo, 'at') as memberships:
    i = 0
    for index, row in P.iterrows():
        if i >= 30: break # only max. 30 pathways
        memberships.write(row['ID'] + "\t" + "\t".join(row['genes']) + "\n")
        i = i + 1
