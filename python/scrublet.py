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

def score_doublets(mtx, doublet_rate=0.1):
    scrub = scr.Scrublet(mtx.T, expected_doublet_rate=doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    return doublet_scores

def call_doublets(doublet_scores, doublet_rate):
    def manual_threshold(scores):
        top_n = int(doublet_rate * len(scores))
        sorted_scores = np.sort(scores)
        threshold = sorted_scores[len(scores)-top_n:].min()
        return threshold
    # scrublet can be conservative -- making sure I get most doublets
    if scrub.threshold_ > 0.3 and sum(predicted_doublets) < (doublet_rate * len(doublet_scores))/2:
        threshold = manual_threshold(doublet_scores)
        if threshold < 0.2:
            threshold = 0.2
        doublets = scrub.call_doublets(threshold=threshold)
    else:
        doublets = scrub.call_doublets()
    return doublets

