import pandas as pd
import os

folder = '/Users/nathandmaulding/Desktop/DREAMIT2/Data/Paul_withregs_v2/diffexp/'

fs = []
tars = []
for each in open('/Users/nathandmaulding/Desktop/DREAMIT/mouse_regulogs.tsv', 'r'):
    e = each.split()
    if e[0] not in fs:
        fs.append(e[0])
    if e[1] not in tars:
        tars.append(e[1])

for file in os.listdir(folder):
    if 'diffexp' not in file:
        pass
    else:
        df = pd.read_csv(folder + file, header=2, index_col=0, sep='\t')
        sig = df[df['FDR'] <= 0.05]
        fsig = sig[sig.index.isin(fs)]
        factors = fsig.shape[0]
        tsig = sig[sig.index.isin(tars)]
        targets = tsig.shape[0]
        branch = file.split('_diffexp.txt')[0]
        print("Branch: ", branch, "Sig Factors: ", factors, "Sig Targets: ", targets)

