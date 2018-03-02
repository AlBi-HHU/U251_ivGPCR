#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os

parser = argparse.ArgumentParser(description = 'Creates a list of the edges')
parser.add_argument('-e', type = str, help = 'name of (complete) edge file 9606.protein.links.v10.HGNC.txt', required = True)
parser.add_argument('-n', type = str, help = 'proteins.nodes file', required = True)
parser.add_argument('-o', type = str, help = 'output file', required = True)

try:
    args = parser.parse_args()
except IOError:
    parser.error()

proteins = []

with open(args.n, 'rt') as nodes:
    reader = csv.DictReader(nodes, delimiter='\t')
    for row in reader:
        proteins.append(row['Identifier'])


with open(args.o, 'wt') as output:
    writer = csv.writer(output, delimiter = '\t')
    with open(args.e, 'rt') as interactome:
        reader = csv.reader(interactome, delimiter='\t')
        for row in reader:
            if row[0] in proteins and row[1] in proteins and row[0] != row[1]: # no self loops!
                writer.writerow([row[0],row[1]])
