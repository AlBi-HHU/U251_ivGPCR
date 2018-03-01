#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
if not os.path.exists("eXamine"):
    os.makedirs("eXamine")

parser = argparse.ArgumentParser(description = 'Transfers heinz output to eXamineS annotations')
parser.add_argument('-m', type = str, help = 'name of file containing module nodes .res', required = True)
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

with open(args.ol, 'wt') as links_output:
    writer = csv.writer(links_output, delimiter='\t')
    for protein in proteins.keys():
        writer.writerow(['small', protein])

with open(args.on, 'wt') as protein_output:
    header = ['Identifier', 'Score', 'Symbol', 'URL']
    writer = csv.DictWriter(protein_output, delimiter = '\t', fieldnames = header)
    writer.writeheader()
    for protein, annotation in proteins.items():
        entry = {}
        entry['Identifier'] = protein
        entry['Score'] = annotation[0]
        if len(annotation) == 2:
            entry['Symbol'] = annotation[1]
        else:
            entry['Symbol'] = protein
        entry['URL'] = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene='+ protein
        writer.writerow(entry)

with open(args.om, 'wt') as annotations:
    header = ['Identifier', 'Category', 'Score', 'Symbol', 'URL']
    line = {'Identifier': 'small', 'Category': 'Module', 'Score': '1', 'Symbol': 'Small', 'URL': 'about:blank'}
    writer = csv.DictWriter(annotations, delimiter = '\t', fieldnames = header)
    writer.writeheader()
    writer.writerow(line)
