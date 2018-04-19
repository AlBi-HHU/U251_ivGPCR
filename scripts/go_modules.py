#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os


parser = argparse.ArgumentParser(description = 'Transfers GO map to eXamineS annotations')
parser.add_argument('-g', type = str, help = 'name of file with GO notations GO_biomart.txt', required=True)
parser.add_argument('-s', type = str, help = 'name of file with GO sets _sets.exm', required=True)
parser.add_argument('-n', type = str, help = 'proteins.nodes file', required = True)
parser.add_argument('-oa', type = str, help = 'go_and_kegg.annotations file', required = True)
parser.add_argument('-ol', type = str, help = 'go_and_kegg.links file', required = True)


try:
    args = parser.parse_args()
except IOError:
    parser.error()

annotations = {}
go_modules = {}
proteins = []

with open(args.n, 'rt') as nodes:
    reader = csv.DictReader(nodes, delimiter='\t')
    for row in reader:
        proteins.append(row['Identifier'])

with open(args.s, 'rt') as sets:
    reader = csv.DictReader(sets, delimiter='\t')
    for row in reader:
        annotations[row['ID']] = ['', row['P-value'], row['Label'], '']
        go_modules[row['ID']] = []

with open(args.g, 'rt') as go_file:
    reader = csv.reader(go_file, delimiter='\t')
    for row in reader:
        if len(row) > 4:
            protein = row[0]
            module = row[1]
            category = row[3]
            if protein in proteins and module in go_modules.keys() and protein not in go_modules[module]:
                go_modules[module].append(protein)
            if module in annotations.keys():
                if category == 'cellular_component':
                    annotations[module][0] = 'Component'
                elif category == 'molecular_function':
                    annotations[module][0] = 'Function'
                elif category == 'biological_process':
                    annotations[module][0] = 'Process'



with open(args.oa, 'wt') as ann_output:
    header = ['Identifier', 'Category', 'Score : Number', 'Symbol', 'URL : Href']
    writer = csv.DictWriter(ann_output, delimiter='\t', fieldnames = header)
    entry = {}
    writer.writeheader()
    for module, annotation in annotations.items():
        if len(go_modules[module]) > 0:
            entry['Identifier'] = module
            entry['Category'] = annotation[0]
            entry['Score : Number'] = annotation[1].replace('< ', '')
            entry['Symbol'] = annotation[2]
            entry['URL : Href'] = 'http://amigo.geneontology.org/cgi-bin/amigo/term_details?term='+ module
            writer.writerow(entry)

with open(args.ol, 'wt') as output:
    writer = csv.writer(output, delimiter = '\t')
    for module in go_modules.keys():
        if len(go_modules[module]) > 0:
            line = [module]
            for el in go_modules[module]:
                line.append(el)
            writer.writerow(line)
