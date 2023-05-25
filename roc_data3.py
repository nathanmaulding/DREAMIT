import pandas as pd
import numpy as np
import os

#folder = '/Users/nathandmaulding/Desktop/1_2/'
#branch = '1_2'
folder = '/projects/sysbio/users/nathanmaulding/results/'

method = 'Pearson'

final_dreamit = {}
final_de = {}
for each in os.listdir(folder):  # traj dirs
    if '_v2' in each:
        group = each
        for b in os.listdir(folder + each):  # branch dirs
            branch = b
            if '.txt' not in b and '_figs' not in b:
                final_dreamit[group + '_' + branch] = {}
                final_de[group + '_' + branch] = {}
                df = pd.DataFrame()
                for f in os.listdir(folder + each + '/' + b):  # files in branch dir
                    if 'dictionary' not in f and branch not in f:
                        if 'Rolling' in method:
                            if method in f:
                                df2 = pd.read_csv(folder + each + '/' + b + '/' + f, header=2, index_col=0, sep='\t')
                                df = pd.concat([df, df2])
                        else:
                            if method in f and 'Rolling' not in f:
                                df2 = pd.read_csv(folder + each + '/' + b + '/' + f, header=2, index_col=0, sep='\t')
                                df = pd.concat([df, df2])
                if df.empty:
                    print('Empty')
                    continue
                df = df.sort_values('FDR')
                df2 = df[~df.index.duplicated(keep='first')]
                df2 = df2.dropna()
                dreamit_branch_fdr = df2['FDR'].tolist()
                dreamit_branch_factor = df2.index.tolist()
                for i in range(len(dreamit_branch_fdr)):
                    final_dreamit[group + '_' + branch][dreamit_branch_factor[i]] = dreamit_branch_fdr[i]
                #final_dreamit = final_dreamit + dreamit_branch_fdr
                de = pd.read_csv(folder + each + '/' + b + '/' + branch + '_diffexp.txt', header=2, index_col=0, sep='\t')
                de2 = de[de.index.isin(df2.index)]
                de_branch_fdr = de2['FDR'].tolist()
                de_branch_factor = de2.index.tolist()
                for j in range(len(de_branch_fdr)):
                    final_de[group + '_' + branch][de_branch_factor[j]] = de_branch_fdr[j]
                #final_de = final_de + de_branch_fdr

dreamit_file = pd.DataFrame(final_dreamit)
dreamit_file.to_csv('/projects/sysbio/users/nathanmaulding/results/dreamit_fdr_results_' + method + '.csv')
de_file = pd.DataFrame(final_de)
de_file.to_csv('/projects/sysbio/users/nathanmaulding/results/de_fdr_results_' + method + '.csv')