#Load required modules
import argparse
import os
import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as smp

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

df = pd.read_csv('KEGG/KEGG_pathways.csv')

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

pvals = list(P['p'])
rej, pcorr = smp.multipletests(pvals, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)[:2]
#print(pvals, pcorr)
P['p_corr'] = pcorr

P = P.sort_values(by=['p_corr'])

with open (args.ao, 'wt') as annotation:
    i = 0
    for index, row in P.iterrows():
        if i >= 30 or row['p_corr'] >= .05: break # only max. 30 pathways
        annotation.write(row['ID'] + "\tKEGG Pathway\t" + str('%.2E' % row['p_corr']) + "\t" + row['description'] + "\thttp://www.kegg.jp\n")
        i = i + 1
        
with open (args.mo, 'wt') as memberships:
    i = 0
    for index, row in P.iterrows():
        if i >= 30 or row['p_corr'] >= .05: break # only max. 30 pathways
        memberships.write(row['ID'] + "\t" + "\t".join(row['genes']) + "\n")
        i = i + 1
