import pandas as pd
import numpy as np
import os

folder = '/projects/sysbio/users/nathanmaulding/results/'

final_genie = {}
for each in os.listdir(folder):
    group = each.split('.csv')[0]
    branch = pd.read_csv(folder + each, index_col=0)
    final_genie[group] = {}
    branch_pval = branch['pvalue'].tolist()
    branch_factor = branch.index.tolist()
    for i in range(len(branch_factor)):
        final_genie[group][branch_factor[i]] = branch_pval[i]

genie_file = pd.DataFrame(final_genie)
genie_file.to_csv(folder + 'genie_fdr_results.csv')