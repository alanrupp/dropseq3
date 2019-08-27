import pandas as pd
import numpy as np
import scrublet as scr
from sklearn.linear_model import LinearRegression

def find_expected_doublet_rate(mtx):
    doublet_rate = pd.read_csv("/home/alanrupp/Programs/dropseq3/data/10x_doublet_rate.csv")
    lm = LinearRegression().fit(np.array(doublet_rate['cells']).reshape(-1,1),
                                np.array(doublet_rate['rate']).reshape(-1,1))
    n_cells = mtx.shape[1]
    expected_doublet_rate = lm.predict(n_cells)[0,0] / 100
    return expected_doublet_rate

def score_doublets(mtx, doublet_rate):
    scrub = scr.Scrublet(mtx.T, expected_doublet_rate=doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    doublets = scrub.call_doublets()
    # if scrublet can't automatically call doublets, take top `doublet_rate` cells
    def manual_doublets(scores):
        top_n = int(doublet_rate * len(scores))
        sorted_scores = np.sort(scores)
        threshold = sorted_scores[len(scores)-top_n:].min()
        doublets = scores > threshold
        return doublets
    if doublets is None:
        doublets = manual_doublets(doublet_scores)
    elif doublet_rate * len(doublets) < sum(doublets):
        doublets = manual_doublets(doublet_scores)
    return doublets
