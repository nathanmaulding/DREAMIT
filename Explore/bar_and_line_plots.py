#####################################################

import numpy as np
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
import pandas as pd



#####################################################
# total, dreamit, de, branchname
f = pd.read_excel('/Users/nathandmaulding/Desktop/DREAMIT_results_table_new.xlsx', header=1)
print(f.columns)

total = f['Database TF Markers'].to_list()
print(total)
dreamit = f['DREAMIT2 TF Markers Unique (Fischer’s or hyper geometric)'].to_list()
print(dreamit)
de = f['DiffExp TF Markers'].to_list()
print(de)
b = f['Dataset'] + '-' + f['Branch']
branchname = b.to_list()
print(branchname)
#fdr001 = f['Unique Significant Pathways (FDR<0.01)'].to_list()
#dc = f['Direct Correlation'].to_list()




## unnormalized ##
plt.bar(branchname, total, color='black', alpha=0.5, label='Total Markers')
plt.bar(branchname, dreamit, color='red', alpha=0.5, label='DREAMIT2')
plt.bar(branchname, de, color='blue', alpha=0.5, label='DE')
plt.xticks(rotation=90)
plt.title('Factor Markers')
plt.legend(loc='upper right')
plt.ylabel('Markers')
plt.xlabel('Dataset')
plt.tight_layout()
plt.savefig('unnormalized_factormarkers.png')
plt.clf()


## normalized ##
ntotal = []
ndreamit = []
nde = []
for i in range(len(total)):
    if total[i] > 0:
        ntotal.append(total[i]/total[i])
        ndreamit.append(dreamit[i]/total[i])
        nde.append(de[i]/total[i])
    else:
        ntotal.append(0)
        ndreamit.append(0)
        nde.append(0)
plt.bar(branchname, ntotal, color='black', alpha=0.5, label='Total Markers')
plt.bar(branchname, ndreamit, color='red', alpha=0.5, label='DREAMIT2')
plt.bar(branchname, nde, color='blue', alpha=0.5, label='DE')
plt.xticks(rotation=90)
plt.title('Factor Markers')
plt.legend(loc='upper right')
plt.ylabel('Marker Proportions')
plt.xlabel('Dataset')
plt.tight_layout()
plt.savefig('normalized_factormarkers.png')
plt.clf()

'''
size = f['Data size'].to_list()
time = f['Runtime (seconds)'].to_list()
zipl = zip(size, time)
sorted_pairs = sorted(zipl)
tuples = zip(*sorted_pairs)
size, time = [list(x) for x in tuples]
plt.plot(size, time, color='red')
plt.scatter(size, time, color='black')
plt.title('Runtime based on data size')
plt.xlabel('Data Size (in MB)')
plt.ylabel('Runtime (in seconds)')
plt.savefig('pathway_runtime.png')
'''
'''
plt.bar(branchname, total, color='black', alpha=0.5, label='Tested')
plt.bar(branchname, dreamit, color='red', alpha=0.5, label='FDR<0.05')
plt.bar(branchname, fdr001, color='green', alpha=0.5, label='FDR<0.01')
plt.bar(branchname, dc, color='yellow', alpha=0.5, label='Direct')
plt.bar(branchname, de, color='blue', alpha=0.5, label='DE')
plt.xticks(rotation=90)
plt.title('Significant Pathways')
plt.legend(loc='upper right', prop={'size': 8})
plt.ylabel('Pathways')
plt.xlabel('Dataset')
plt.tight_layout()
plt.savefig('pathways.png')
plt.clf()
'''
