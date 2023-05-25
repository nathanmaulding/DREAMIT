import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.stats import gaussian_kde
from matplotlib.gridspec import GridSpec
import json


def openFile(path, ivar):
    data = open(path, 'r')
    expDict = ''
    for i, line in enumerate(data):
        line.split()
        if i == ivar:
            expDict = line
    expDict = eval(expDict)
    return expDict

path2json = '/Users/nathandmaulding/Desktop/DREAMIT/Data/Paul_withregs_v2/lucas_paul15.json'
path2 = '/Users/nathandmaulding/Desktop/v20Data/paul_withregs/stem_to_erythro_Smoothed_Expression_dictionary.txt'
branch = 'stem_to_erythro'
b = '_stem_to_erythro'
gene_list2 = ['E2f4', 'Myc']
dataset = 'paul_withregs'
var = 'psuedotime'
exp2 = openFile(path2, 3)

with open(path2json) as j:
    exp = json.load(j)


rawexp = pd.DataFrame(exp['branches'][branch])
smoothexp = pd.DataFrame(exp2)
print(rawexp.columns.tolist())
print(exp.keys())
gene_list = []
for gene in gene_list2:
    if gene in rawexp.columns.tolist():
        gene_list.append(gene)
for gene in gene_list:
    rawptime = [float(i) for i in rawexp[var].tolist()]
    generaw = [float(i) if float(i) > 0 else 0.000000001 for i in rawexp[gene].tolist()]
    # smoothptime needs to be changed according to bin number
    print(gene, len(smoothexp[gene].tolist()), smoothexp[gene].tolist())
    print(gene, len(generaw), len(rawptime), generaw)
    print(rawptime)
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
    #fig = plt.figure()
    # Calculate the point density
    xy = np.vstack([rawptime, generaw])
    z = gaussian_kde(xy)(xy)

    fig = plt.figure()
    gs = GridSpec(4,4)
    ax_joint = fig.add_subplot(gs[1:4, 0:3])
    ax_marg_x = fig.add_subplot(gs[0, 0:3])
    ax_marg_y = fig.add_subplot(gs[1:4, 3])
    uplim = min(smoothexp['psuedotime_med'].tolist()) - min(rawptime)
    ax_joint.set_xlim([min(rawptime), max(smoothexp['psuedotime_med'].tolist())+uplim])
    ax_marg_x.set_xlim([min(rawptime), max(smoothexp['psuedotime_med'].tolist())+uplim])
    ax_joint.scatter(rawptime, generaw, c=z, s=8)
    ax_joint.plot(smoothexp['psuedotime_med'].tolist(), smoothexp[gene].tolist(), color='black', linewidth=2, alpha=0.65)
    ax_joint.scatter(smoothexp['psuedotime_med'].tolist(), smoothexp[gene].tolist(), color='black', s=25)
    ##### set range = xlim setting if used for below #####
    ax_marg_x.hist(rawptime, range=[min(rawptime), max(smoothexp['psuedotime_med'].tolist())+uplim], bins=20, color='gray')
    #######################################################
    ax_marg_y.hist(generaw, orientation='horizontal', bins=20, color='gray')
    # Set labels on joint
    ax_joint.set_xlabel('Pseudotime')
    ax_joint.set_ylabel('Normalized Expression')

    # Turn off tick labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)

    # Set labels on marginals
    ax_marg_y.set_xlabel('No. Cells')
    ax_marg_x.set_ylabel('No. Cells')
    #plt.title(gene + ' ' + branch + ' raw and smooth expression')
    plt.savefig('./figs/' + gene + '_' + dataset + b + '.png', dpi=300)
    plt.show()
    print(gene)
    print(smoothexp[gene].tolist())
    print(np.log(smoothexp[gene].tolist()))
    print(min(generaw))
    print(np.log(min(generaw)))

    plt.plot(smoothexp['psuedotime_med'].tolist(), smoothexp[gene].tolist(), color='black', linewidth=2, alpha=0.65)
    #plt.xlim([0, 200])
    plt.xlabel('Pseudotime')
    plt.ylabel('Expression')
    plt.savefig('./figs/' + gene + '_' + dataset + b + '_smoothexp.png', dpi=300)
    plt.show()

average_exp = []
for gene in rawexp.columns.tolist()[9:]:
    #print(rawexp[gene].tolist())
    average_exp.append(np.mean([float(i) for i in rawexp[gene].tolist()]))
print(average_exp)
print(np.mean(average_exp))
p = sorted([float(i) for i in rawexp[var].tolist()])
Q1, Q3 = np.percentile(p, [25, 75])
IQR = Q3 - Q1
ul = Q3 + 1.5 * IQR
ll = Q1 - 1.5 * IQR
print(ll, ul)

rawexp[var] = [float(i) for i in rawexp[var].tolist()]
rawexp = rawexp[(rawexp[var] < ul) & (rawexp[var] > ll)]
print(len(rawexp[var]), rawexp)




'''
    #plt.scatter(rawptime, np.log(generaw), c=z, s=8)
    #plt.plot(smoothexp['psuedotime_med'].tolist(), np.log(smoothexp[gene].tolist()), color='black', linewidth=2, alpha=0.65)
    #plt.scatter(smoothexp['psuedotime_med'].tolist(), np.log(smoothexp[gene].tolist()), color='black', s=25)
    print(gene)
    print(smoothexp[gene].tolist())
    print(np.log(smoothexp[gene].tolist()))
    #plt.plot(x2, splineexp, color='red', linewidth=4)
    plt.xlabel('pseudotime')
    plt.ylabel('log expression')
    plt.xlim([0,200])
    #plt.xticks(np.arange(min(rawptime), max(rawptime)+((max(rawptime)-min(rawptime))/10), step=(max(rawptime)-min(rawptime))/10))
    plt.title(gene + ' ' + branch + ' raw expression')
    plt.savefig(gene + '_' + dataset + '.png', dpi=300)
    plt.show()
'''
