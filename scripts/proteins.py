#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
if not os.path.exists("eXamine"):
    os.makedirs("eXamine")

parser = argparse.ArgumentParser(description = 'Transfers heinz output to eXamine annotations')
parser.add_argument('-m', type = str, help = 'name of file containing module nodes .res', required = True)
parser.add_argument('-fc', type = str, help = 'name of file containing fold change values', required = True)
parser.add_argument('-ol', type = str, help = 'links file name', required = True)
parser.add_argument('-on', type = str, help = 'nodes file name', required = True)
parser.add_argument('-om', type = str, help = 'modules annotation file name', required = True)

try:
    args = parser.parse_args()
except IOError:
    parser.error()

proteins = {}

with open(args.m, 'rt') as heinzfile:
    reader = csv.DictReader(heinzfile, delimiter='\t',fieldnames = ( "#label","score"))
    for row in reader:
        proteins[row['#label']] = [row['score']]

#print(proteins)

with open(args.fc, 'rt') as foldchangefile:
    reader = csv.DictReader(foldchangefile, delimiter='\t')#, fieldnames = ("Id", "log2FoldChange", "padj"), header=0)
    for row in reader:
        proteins[row['']].append(row['log2FoldChange'])
        #proteins[row['']].append(row['padj'])
        #print(proteins[row['Id']])#.append(row['log2FoldChanage'])

#print(proteins)

with open(args.ol, 'wt') as links_output:
    writer = csv.writer(links_output, delimiter='\t')
    for protein in proteins.keys():
        writer.writerow(['small', protein])

with open(args.on, 'wt') as protein_output:
    header = ['Identifier', 'FC : Number', 'Score : Number', 'Symbol', 'URL : Href']
    writer = csv.DictWriter(protein_output, delimiter = '\t', fieldnames = header)
    writer.writeheader()
    for protein, annotation in proteins.items():
        #print(protein, annotation)
        entry = {}
        entry['Identifier'] = protein
        entry['FC : Number'] = annotation[1]
        entry['Score : Number'] = annotation[0]
        if len(annotation) == 3:
            entry['Symbol'] = annotation[3]
        else:
            entry['Symbol'] = protein
        entry['URL : Href'] = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene='+ protein
        writer.writerow(entry)

with open(args.om, 'wt') as annotations:
    header = ['Identifier', 'Category', 'Score', 'Symbol', 'URL']
    line = {'Identifier': 'small', 'Category': 'Module', 'Score': '1', 'Symbol': 'Small', 'URL': 'about:blank'}
    writer = csv.DictWriter(annotations, delimiter = '\t', fieldnames = header)
    writer.writeheader()
    writer.writerow(line)
