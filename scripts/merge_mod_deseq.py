#!/usr/bin/env python

#Load required modules
import sys, os, argparse
import pandas as pd

# Parse arguments
def get_parser():
    description = 'Merge heinz results with deseq2'
    parser = argparse.ArgumentParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-m', '--module', required=True, help='File name of Heinz module')
    parser.add_argument('-d', '--deseq2', required=True, help='File name of deseq2 results')
    parser.add_argument('-o', '--output', required=True, help='File name for merged results')
    return parser

if __name__ == '__main__':
    # Parse arguments
    args = get_parser().parse_args( sys.argv[1:] )

    try:
        module = pd.read_csv(args.module, sep='\t', header=-1, index_col=0)#[0]
        #print(module)
        deseq2 = pd.read_csv(args.deseq2, header=0, sep='\t', index_col=0)#[['log2FoldChange', 'padj']] # get gene name, log2FC, padj
        #print(deseq2)
        #print(module.columns.values)
        #print(deseq2.columns.values)

        merged_results = pd.merge(module, deseq2, left_index=True, right_index=True, sort=False, how='inner')
        #print(merged_results[['log2FoldChange', 'padj']])
        # win_rev = sum(merged_results['Final IC'] < merged_results['IC'].round(2))
        # win_sup = sum(merged_results['Final IC'] > merged_results['IC'].round(2))
        # draw = sum(merged_results['Final IC'] == merged_results['IC'].round(2))
        # print win_rev, draw, win_sup
        #
        #merged_results.rename(columns = {'1': 'sds'}, inplace=True)
        merged_results[['log2FoldChange', 'padj']].to_csv(args.output, sep='\t')
    except:
        print("merging heinz with deseq2 failed")
