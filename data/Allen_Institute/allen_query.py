import os
import requests
import pandas as pd
import numpy as np
import argparse

# grab command line input
parser = argparse.ArgumentParser()
parser.add_argument('--queries', nargs='*', help='Enter queries separated by spaces')
parser.add_argument('--contrast', help='Enter a single contrast structure')
parser.add_argument('--coronal', action='store_true', help='Only coronal (default: False)')
parser.add_argument('--threshold', nargs=1, type=int, default=1, \
                    help='Expression threshold for query')
args = parser.parse_args()

# filter area information by matches for query and contrast
def filter_areas(areas, contrast):
    all_areas = pd.read_csv('structures.tsv', sep='\t')
    hit = list()
    for i, df in all_areas.iterrows():
        match_query = len(set(df.values).intersection(set(areas)))
        match_contrast = len(set(df.values).intersection(set([contrast])))
        if match_query + match_contrast == 2:
            hit.append(True)
        else:
            hit.append(False)
    return(all_areas[hit])

# find stucture ID by contrast
def find_contrast(contrast):
    all_areas = pd.read_csv('structures.tsv', sep='\t')
    # keep only rows that list contrast
    contrast_match = list()
    for i, df in all_areas.iterrows():
        if len(set(df.values).intersection(set([contrast]))) > 0:
            contrast_match.append(True)
        else:
            contrast_match.append(False)
    all_areas = all_areas[contrast_match]
    all_areas = all_areas.reset_index()
    # keep row with most np.nan
    nans = list()
    for i, df in all_areas.iterrows():
        nans.append(int(sum(pd.isnull(df.values))))
    # find max position in list
    max_pos = nans.index(max(nans))
    return(all_areas.iloc[max_pos]['structureID'])

# generate request text
def make_request(query, contrast, threshold):
    request = "http://api.brain-map.org/api/v2/data/query.csv?criteria=service::mouse_differential"
    if args.coronal:
        request += "[set$eq'mouse_coronal']"
    else:
        request += "[set$eq'mouse']"
    request += "[structures1$eq{}]".format(contrast)
    request += "[structures2$eq{}]".format(query)
    request += "[threshold2$eq{},50]".format(threshold)
    return(request)

# call request
def get_request(request):
    return(requests.get(request).text)

# write data to CSV
def write_csv(request, fname):
    with open(fname, 'w') as f:
        f.write(request)

# - Run program ---------------------------------------------------------------
if __name__ == '__main__':
    queries = filter_areas(args.queries, args.contrast)
    contrast_id = find_contrast(args.contrast)
    contrast_name = args.contrast
    print(args)

    for i, df in queries.iterrows():
        print('grabbing ' + df.iloc[-1] + ' vs. ' + contrast_name)
        r = make_request(query = df['structureID'], contrast = contrast_id, threshold = args.threshold)
        r = get_request(r)
        fname = contrast_name + "/" + df.iloc[-1] + ".csv"
        write_csv(r, fname)
