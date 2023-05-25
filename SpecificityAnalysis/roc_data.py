import pandas as pd
import numpy as np
import os

#folder = '/Users/nathandmaulding/Desktop/1_2/'
#branch = '1_2'
folder = '/projects/sysbio/users/nathanmaulding/results/'

final_dreamit = []
final_de = []
for each in os.listdir(folder):
    if '_v2' in each:
        group = each
        for b in os.listdir(folder + each):
            branch = b
            if '.txt' not in b and '_figs' not in b:
                df = pd.DataFrame()
                for f in os.listdir(folder + each + '/' + b):
                    if 'dictionary' not in f and branch not in f:
                        df2 = pd.read_csv(folder + each + '/' + b + '/' + f, header=2, index_col=0, sep='\t')
                        df = pd.concat([df, df2])
                df = df.sort_values('FDR')
                df2 = df[~df.index.duplicated(keep='first')]
                df2 = df2.dropna()
                dreamit_branch_fdr = df2['FDR'].tolist()
                final_dreamit = final_dreamit + dreamit_branch_fdr
                de = pd.read_csv(folder + each + '/' + b + '/' + branch + '_diffexp.txt', header=2, index_col=0, sep='\t')
                de2 = de[de.index.isin(df2.index)]
                de_branch_fdr = de2['FDR'].tolist()
                final_de = final_de + de_branch_fdr

dreamitfile = open('dreamit_fdr_results.txt', 'w')
for e in final_dreamit:
    e2 = str(e) + '\t'
    dreamitfile.write(e2)
dreamitfile.close()

defile = open('de_fdr_results.txt', 'w')
for i in final_de:
    i2 = str(i) + '\t'
    defile.write(i2)
defile.close()


'''
df = pd.DataFrame()
for f in os.listdir(folder):
    if 'dictionary' not in f and branch not in f:
        df2 = pd.read_csv(folder + f, header=2, index_col=0, sep='\t')
        df = pd.concat([df, df2])
df = df.sort_values('FDR')
df2 = df[~df.index.duplicated(keep='first')]
df2 = df2.dropna()
dreamit_branch_fdr = df2['FDR'].tolist()
de = pd.read_csv(folder + branch + '_diffexp.txt', header=2, index_col=0, sep='\t')
de2 = de[de.index.isin(df2.index)]
de_branch_fdr = de2['FDR'].tolist()
'''
