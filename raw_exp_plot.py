import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate


def openFile(path, ivar):
    data = open(path, 'r')
    expDict = ''
    for i, line in enumerate(data):
        line.split()
        if i == ivar:
            expDict = line
    expDict = eval(expDict)
    return expDict

#path = '/Users/nathandmaulding/Desktop/bin52_ebi10_raw.txt'
#path2 = '/Users/nathandmaulding/Desktop/ExpresionPlots/bin52/1_2_Smoothed_Expression_dictionary.txt'
path = '/Users/nathandmaulding/Desktop/test_binopt2_ebi10_raw.txt'
path2 = '/Users/nathandmaulding/Desktop/test_binopt2_Smoothed_Expression_dictionary.txt'
branch = '1 -> 2'
gene_list = ['MEF2C', 'WT1', 'ELK1']
dataset = 'test_binopt2'
var = 'psuedotime'
exp = openFile(path, 2)
exp2 = openFile(path2, 3)


rawexp = pd.DataFrame(exp[branch])
smoothexp = pd.DataFrame(exp2)
print(rawexp.columns.tolist())
print(exp.keys())
for gene in gene_list:
    rawptime = [float(i) for i in rawexp[var].tolist()]
    generaw = [float(i) for i in rawexp[gene].tolist()]
    # smoothptime needs to be changed according to bin number
    print(gene, len(smoothexp[gene].tolist()), smoothexp[gene].tolist())
    smoothptime = [i for i in np.arange(min(rawptime), max(rawptime), (max(rawptime)-min(rawptime))/len(smoothexp[gene].tolist()))]
    ############
    x = rawptime
    y = generaw
    pairs = []
    for i in range(len(x)):
        if (x[i], y[i]) not in pairs:
            pairs.append((x[i], y[i]))
    x, y = [], []
    for j in pairs:
        x.append(j[0])
        y.append(j[1])
    x, y = zip(*sorted(zip(x, y)))
    x2 = []
    y2 = []
    y3 = []
    for i in range(len(x)):
        if x[i] not in x2:
            if len(y3) > 0:
                y2.append(np.mean(y3))
                y3 = []
            x2.append(x[i])
            y3.append(y[i])
        else:
            y3.append(y[i])
    y2.append(np.mean(y3))
    spl = interpolate.UnivariateSpline(x2, y2)
    splineexp = spl(np.array(x2))
    ################
    plt.scatter(rawptime, generaw)
    plt.plot(smoothptime, smoothexp[gene].tolist(), color='black', linewidth=4)
    #plt.plot(x2, splineexp, color='red', linewidth=4)
    plt.xlabel('pseudotime')
    plt.ylabel('expression')
    plt.xticks(np.arange(min(rawptime), max(rawptime)+((max(rawptime)-min(rawptime))/10), step=(max(rawptime)-min(rawptime))/10))
    plt.title(gene + ' ' + branch + ' raw expression')
    plt.savefig(gene + '_' + dataset + '.png', dpi=300)
    plt.show()
