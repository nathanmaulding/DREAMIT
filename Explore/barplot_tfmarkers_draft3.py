import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import os


detars_df = pd.read_csv('/Users/nathandmaulding/Desktop/specificity_results/v20/figs_allmarkers/detargets/summary.csv',index_col=0)
de_df = pd.read_csv('/Users/nathandmaulding/Desktop/specificity_results/v20/figs_allmarkers/de/summary.csv',index_col=0)
genie_df = pd.read_csv('/Users/nathandmaulding/Desktop/specificity_results/v20/figs_allmarkers/genie/summary.csv', index_col=0)

study_var = 'F1 AUC'

detars_df = detars_df[[col for col in detars_df if study_var in col]]
de_df = de_df[[col for col in de_df if study_var in col]]
genie_df = genie_df[[col for col in genie_df if study_var in col]]


df = detars_df
df = df.join(de_df[[study_var+' rawDE', study_var+' smoothDE']])
df.index = genie_df.index
df = df.join(genie_df[[study_var+' genie', study_var+' smooth_genie']])

print(df)

dft = df.T

for branch in dft.columns.tolist():
    for i in range(len(dft[branch])):
        print(dft[branch][i], dft.index[i], branch)
        plt.bar(i+1, dft[branch][i])
    plt.ylabel(study_var)
    plt.xticks(range(1, len(dft[branch])+1), dft.index.tolist(), rotation=90)
    plt.title(branch)
    plt.show()
