#!/usr/bin/python

'''
Using the CELLEX framework for finding cell-type specificity of genes
'''

import numpy as np
import pandas as pd
import cellex

def run_cellex(data, metadata):
    eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)
    eso.compute(verbose=True)
    return eso.results['esmu']

if __name__ == '__main__':
    # parse command line arguments
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-data', default='counts.csv', help='counts CSV file')
    parser.add_argument('-metadata', default='clusters.csv', help='metadata CSV file')
    args = parser.parse_args()
    # read in data
    counts = pd.read_csv(args.data, index_col=0)
    clusters = pd.read_csv(args.metadata, index_col=0)
    # run CELLEX and return ESmu
    result = run_cellex(counts, clusters)
