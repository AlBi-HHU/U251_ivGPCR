#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description = 'Removes a module from a network')
parser.add_argument('-n', type = str, help = 'network', required = True)
parser.add_argument('-m', type = str, help = 'module.mod', required = True)
parser.add_argument('-o', type = str, help = 'output network', required = True)

try:
    args = parser.parse_args()
except IOError:
    parser.error()

nodes = {}

try:
    module_table = pd.read_csv(args.m, sep='\t', header=None)
    nodes = set(module_table[0])
    network = pd.read_csv(args.n, sep='\t', header=None)
    for index, e in network.iterrows():
        if e[0] in nodes or e[1] in nodes: network.drop(index, inplace=True)
    network.to_csv(args.o, sep='\t', index=False, header=False)
except:
    print("some error occured")
