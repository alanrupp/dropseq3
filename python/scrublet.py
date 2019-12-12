#!/usr/bin/python

import pandas as pd
import numpy as np
import scrublet as scr
from sklearn.linear_model import LinearRegression

def find_expected_doublet_rate(mtx):
    doublet_rate = pd.read_csv("/home/alanrupp/Programs/dropseq3/data/10x_doublet_rate.csv")
    lm = LinearRegression().fit(np.array(doublet_rate['cells']).reshape(-1,1),
                                np.array(doublet_rate['rate']).reshape(-1,1))
    n_cells = mtx.shape[1]
    expected_doublet_rate = lm.predict(np.array(n_cells).reshape(-1,1))[0,0] / 100
    return expected_doublet_rate

def score_doublets(mtx, doublet_rate):
    scrub = scr.Scrublet(mtx.T, expected_doublet_rate=doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    def manual_threshold(scores):
        top_n = int(doublet_rate * len(scores))
        sorted_scores = np.sort(scores)
        threshold = sorted_scores[len(scores)-top_n:].min()
        return threshold
    # scrublet can be conservative -- making sure I get most doublets
    if scrub.threshold_ > 0.3 and sum(predicted_doublets) < (doublet_rate * len(doublet_scores))/2:
        threshold = manual_threshold(doublet_scores)
        doublets = scrub.call_doublets(threshold=threshold)
    else:
        doublets = scrub.call_doublets()
    return doublets

if __name__ == '__main__':
    from scipy import io
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='matrix folder location')
    args = parser.parse_args()

    # import mtx
    mtx = io.mmread(f"{path}/matrix.mtx.gz")

    # calcualte doublet rate and find doublets
    doublet_rate = find_expected_doublet_rate(mtx)
    doublets = score_doublets(mtx, doublet_rate)
