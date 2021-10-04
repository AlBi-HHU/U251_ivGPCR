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

signal_transduction_32 = ["path:hsa02020", "path:hsa04014", "path:hsa04015", "path:hsa04010", "path:hsa04013", "path:hsa04016", "path:hsa04011", "path:hsa04012", "path:hsa04310", "path:hsa04330", "path:hsa04340", "path:hsa04341", "path:hsa04350", "path:hsa04390", "path:hsa04392", "path:hsa04392", "path:hsa04370", "path:hsa04371", "path:hsa04630", "path:hsa04064", "path:hsa04668", "path:hsa04066", "path:hsa04068", "path:hsa04020", "path:hsa04070", "path:hsa04072", "path:hsa04071", "path:hsa04024", "path:hsa04022", "path:hsa04151", "path:hsa04152", "path:hsa04150", "path:hsa04075"]

signaling_mol_33 = ["path:hsa04080", "path:hsa04060", "path:hsa04512", "path:hsa04514"]

cellular_processes_4 = ["path:hsa04144", "path:hsa04145", "path:hsa04142", "path:hsa04146", "path:hsa04140", "path:hsa04138", "path:hsa04137",
"path:hsa04110", "path:hsa04210", "path:hsa04216", "path:hsa04217", "path:hsa04115", "path:hsa04218", "path:hsa04510", "path:hsa04520", "path:hsa04530", "path:hsa04540", "path:hsa04550"]

diseases_6 = ["path:hsa05200", "path:hsa05230", "path:hsa05231", "path:hsa05202", "path:hsa05206", "path:hsa05205", "path:hsa05204", "path:hsa05203", "path:hsa05210", "path:hsa05212", "path:hsa05225", "path:hsa05226", "path:hsa05214", "path:hsa05216", "path:hsa05221", "path:hsa05220", "path:hsa05217", "path:hsa05218", "path:hsa05211", "path:hsa05219", "path:hsa05215", "path:hsa05213", "path:hsa05224", "path:hsa05222", "path:hsa05223"]

extra = ["path:hsa01100", "path:hsa04550", "path:hsa"]


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
    #print(module, pathway, module_pathway)
    module_non_pathway = len(module - pathway)
    non_module_pathway = len(pathway - module)
    non_module_non_pathway = len(background - (module | pathway))

    #print(module_pathway, module_non_pathway, non_module_pathway, non_module_non_pathway)

    oddsratio, pvalue = stats.fisher_exact([[module_pathway, non_module_pathway], [module_non_pathway, non_module_non_pathway]])
    #oddsratio, pvalue = stats.fisher_exact([[1, 9], [11, 3]])
    #print(pvalue)
    #P.append(pvalue, row['ID'], row['description'], list(module & pathway))
    P.at[index, 'p'] = pvalue
    P.at[index, 'shared'] = module_pathway
    P.at[index, 'ID'] = row['ID']
    P.at[index, 'description'] = row['description'][:-23]
    P.at[index, 'genes'] = list(module & pathway)

    #print(pvalue, module_pathway, module_non_pathway, non_module_pathway, non_module_non_pathway, row['description'])

pvals = list(P['p'])
rej, pcorr = smp.multipletests(pvals, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)[:2]
#print(pvals, pcorr)

#P.to_csv("foo.csv")

P['p_corr'] = pcorr

P = P.sort_values(by=['p_corr'])
#print(P)

with open (args.ao, 'at') as annotation:
    #add enriched pathways
    i = 0
    for index, row in P.iterrows():
        #print(row)
        if i >= 30 or row['p_corr'] >= .05: break # only max. 30 pathways
        annotation.write("KEGG_all_" + row['ID'] + "\tKEGG all enriched\t" + str('%.2E' % row['p_corr']) + "\t" + row['description'] + "\thttp://www.kegg.jp\n")
        i = i + 1
    annotation.flush()
    #print("now all")
    # add signaling pathways
    for index, row in P.iterrows():
        #print(row)
        if row['ID'] in signal_transduction_32 and row['shared'] > 0:
            annotation.write("KEGG_ST_" + row['ID'] + "\tKEGG signal transduction\t" + str('%.2E' % row['p_corr']) + "\t" + row['description'] + "\thttp://www.kegg.jp\n")
        annotation.flush()
        if row['ID'] in signaling_mol_33 and row['shared'] > 0:
            annotation.write("KEGG_SMI_" + row['ID'] + "\tKEGG signaling mol. and int.\t" + str('%.2E' % row['p_corr']) + "\t" + row['description'] + "\thttp://www.kegg.jp\n")
        annotation.flush()
        if row['ID'] in cellular_processes_4 and row['shared'] > 0:
            annotation.write("KEGG_CP_" + row['ID'] + "\tKEGG cellular processes\t" + str('%.2E' % row['p_corr']) + "\t" + row['description'] + "\thttp://www.kegg.jp\n")
        annotation.flush()
        if row['ID'] in diseases_6 and row['shared'] > 0:
            annotation.write("KEGG_dis_" + row['ID'] + "\tKEGG diseases\t" + str('%.2E' % row['p_corr']) + "\t" + row['description'] + "\thttp://www.kegg.jp\n")
        annotation.flush()
        if row['ID'] in extra and row['shared'] > 0:
            annotation.write("KEGG_extra_" + row['ID'] + "\tKEGG extra\t" + str('%.2E' % row['p_corr']) + "\t" + row['description'] + "\thttp://www.kegg.jp\n")



with open (args.mo, 'at') as memberships:
    #i = 0
    for index, row in P.iterrows():
        #if i >= 30 or row['p_corr'] >= .05: break # only max. 30 pathways
        memberships.write("KEGG_all_" + row['ID'] + "\t" + "\t".join(row['genes']) + "\n")
        memberships.write("KEGG_ST_" + row['ID'] + "\t" + "\t".join(row['genes']) + "\n")
        memberships.write("KEGG_SMI_" + row['ID'] + "\t" + "\t".join(row['genes']) + "\n")
        memberships.write("KEGG_CP_" + row['ID'] + "\t" + "\t".join(row['genes']) + "\n")
        memberships.write("KEGG_dis_" + row['ID'] + "\t" + "\t".join(row['genes']) + "\n")
        memberships.write("KEGG_extra_" + row['ID'] + "\t" + "\t".join(row['genes']) + "\n")
        #i = i + 1
