#!/usr/bin/python

'''
Using the CELLEX framework for finding cell-type specificity of genes
'''

import numpy as np
import pandas as pd
import cellex

def run_cellex(data, metadata, human=False, result='esmu'):
    eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)
    eso.compute(verbose=True)
    if human:
        cellex.utils.mapping.mouse_symbol_to_mouse_ens(eso.results["esmu"], drop_unmapped=True, verbose=True)
        cellex.utils.mapping.mouse_ens_to_human_ens(eso.results["esmu"], drop_unmapped=True, verbose=True)
    return eso.results[result]
