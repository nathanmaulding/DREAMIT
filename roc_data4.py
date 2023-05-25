import pandas as pd
import numpy as np
import os

folder = '/projects/sysbio/users/nathanmaulding/results/'
ks_thresh = 0.5

final_dreamit = {}
final_de = {}
for each in os.listdir(folder):
    if '_v2' in each:
        group = each
        for b in os.listdir(folder + each):
            branch = b
            if '.txt' not in b and '_figs' not in b:
                final_dreamit[group + '_' + branch] = {}
                final_de[group + '_' + branch] = {}
                df = pd.DataFrame()
                for f in os.listdir(folder + each + '/' + b):
                    if 'dictionary' not in f and branch not in f:
                        df2 = pd.read_csv(folder + each + '/' + b + '/' + f, header=2, index_col=0, sep='\t')
                        df = pd.concat([df, df2])
                df = df.sort_values('pvalue')
                print('before', df)
                df = df[df['KS Statistic Value'] >= ks_thresh]
                print('after', df)
                df2 = df[~df.index.duplicated(keep='first')]
                df2 = df2.dropna()
                dreamit_branch_fdr = df2['pvalue'].tolist()
                dreamit_branch_factor = df2.index.tolist()
                for i in range(len(dreamit_branch_fdr)):
                    final_dreamit[group + '_' + branch][dreamit_branch_factor[i]] = dreamit_branch_fdr[i]

dreamit_file = pd.DataFrame(final_dreamit)
dreamit_file.to_csv('/projects/sysbio/users/nathanmaulding/results/dreamit_fdr_results4.csv')