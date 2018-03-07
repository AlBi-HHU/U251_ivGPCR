#Load required modules
import sys, os

from Bio.KEGG import REST

import pandas as pd

human_pathways = REST.kegg_list("pathway", "hsa").read()

genes = {}
description = {}

df = pd.DataFrame(columns = ['ID', 'description', 'genes'])

for i, p in enumerate(human_pathways.rstrip().split("\n")):
    p_entry, p_description = p.split("\t")
    genes[p_entry] = []
    df.at[i, 'ID'] = p_entry
    df.at[i, 'description'] = p_description

    description[p_entry] = p_description
    #print(p_entry, p_description)
    pathway_file = REST.kegg_get(p_entry).read()
    # iterate through each KEGG pathway file, keeping track of which section
    # of the file we're in, only read the gene in each pathway
    current_section = None
    genelist = []
    for line in pathway_file.rstrip().split("\n"):

        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section

        if current_section == "GENE":
            try:
                gene_identifiers, gene_description = line[12:].split("; ")
                gene_id, gene_symbol = gene_identifiers.split()

                genelist.append(gene_symbol)
            except:
                print("WARNING -- pathway %s missing a gene with unknown gene symbol" % p_entry)
    df.at[i, 'genes'] = genelist
    #print(p_entry, genes[p_entry])

df.to_csv('KEGG_pathways.csv')
    
    

