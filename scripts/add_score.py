#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description = 'Adds a score to a .nodes file')
parser.add_argument('-n', type = str, help = 'nodes file name', required = True)
parser.add_argument('-s', type = str, help = 'name of file containing additional score', required = True)

try:
    args = parser.parse_args()
except IOError:
    parser.error()

try:
    nodes_table = pd.read_csv(args.n, sep='\t', header=0, index_col=0)
    #print(nodes_table)
    add_score_table = pd.read_csv(args.s, sep='\t', header=0, index_col=0)
    add_score_table.drop(add_score_table.columns[[0, 2, 3, 4, 5]], axis=1, inplace=True)  # df.columns is zero-based pd.Index
    add_score_table.rename(columns={'log2FoldChange': 'FC other : Number'}, inplace=True)
    merged_table = pd.merge(nodes_table, add_score_table, left_index=True, right_index=True, sort=False, how='inner')
    merged_table.index.name = 'Identifier'
    #print(merged_table)
    merged_table.to_csv(args.n, sep='\t')
except:
    print("adding score failed")
