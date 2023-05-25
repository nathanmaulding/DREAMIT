# this script will plot the gene expression of the optimal param

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

aicdf = pd.read_csv('/Users/nathandmaulding/Desktop/splinechoice/EBI5_2_3_aic_df.csv', index_col=0) #pd.read_csv('/Users/nathandmaulding/Desktop/ebi5_3_4_aic_df.csv', index_col=0)
print(aicdf)
aiccols = []
for col in aicdf.columns:
    aiccols.append(int(col))
aicdf.columns = aiccols

s = open('/Users/nathandmaulding/Desktop/splinechoice/EBI5_2_3_subsample_dictionary.txt', 'r') #open('/Users/nathandmaulding/Desktop/small_sub2.txt', 'r')
sub = ''
for i, line in enumerate(s):
    if i > 0:
        sub = sub + line
subdict = eval(sub)
print(subdict.keys())

# Best - Branch 2 -> 3 Smoothing Factor of  7.5 Bin number of  5 - CV of 0.81
# Worst - Branch 2 -> 3 Smoothing Factor of  0.10 Bin number of  29 - CV of 22.92
genelist = ['FRY']
title = 'EBI5_2_3_'

for g in genelist:
    gene = subdict[g]

    bins = 29
    sf = 0.10
    bin_inc = (1 / bins) / 2
    for exp in gene[bins][sf][3]:  # [bin][sf][exp index is 3]
        ptime = [(i/bins)-bin_inc for i in range(1, bins+1)]
        plt.plot(ptime, exp, alpha=0.5)
    plt.ylabel('Expression')
    plt.xlabel('Psuedotime')
    plt.title(g + ' with worst params')
    plt.savefig(title + g + '_worst_params.png', dpi=300)
    plt.show()

    bins = 5
    sf = 7.5
    bin_inc = (1 / bins) / 2
    for exp in gene[bins][sf][3]:
        ptime = [(i/bins)-bin_inc for i in range(1, bins+1)]
        plt.plot(ptime, exp, alpha=0.5)
    plt.ylabel('Expression')
    plt.xlabel('Psuedotime')
    plt.title(g + ' with best params')
    plt.savefig(title + g + '_best_params.png', dpi=300)
    plt.show()
