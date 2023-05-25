import pandas as pd
import random
import numpy as np
from scipy.stats import ks_2samp
import scipy.stats as stats
import os

genie_dir = '/projects/sysbio/users/nathanmaulding/results/genie3/'

f = open('/projects/sysbio/users/nathanmaulding/data/trrust_rawdata.human.tsv', 'r')
#f = open('/Users/nathandmaulding/Desktop/DREAMIT2/trrust_rawdata.human.tsv', 'r')
tfs = {}
added = {}
doubles = {}
direction_dict = {'Activation': '+', 'Repression': '-', 'Unknown': '?'}
for eachline in f:
    each = eachline.split()
    if each[0] in added.keys():
        if each[1] in added[each[0]]:
            doubles[each[0]].append(each[1])
        else:
            added[each[0]].append(each[1])
    else:
        added[each[0]] = []
        added[each[0]].append(each[1])
        doubles[each[0]] = []

f2 = open('/projects/sysbio/users/nathanmaulding/data/trrust_rawdata.human.tsv', 'r')
#f2 = open('/Users/nathandmaulding/Desktop/DREAMIT2/trrust_rawdata.human.tsv', 'r')
for eachline2 in f2:
    each = eachline2.split()
    if each[0] in tfs.keys():
        if each[0] != each[1]:  # remove auto-correlation
            if each[1] in doubles[each[0]]:
                if [each[1], '?'] not in tfs[each[0]]:
                    tfs[each[0]].append([each[1], '?'])
            else:
                tfs[each[0]].append([each[1], direction_dict[each[2]]])
    else:
        tfs[each[0]] = []
        if each[0] != each[1]:  # remove auto-correlation
            if each[1] in doubles[each[0]]:
                if [each[1], '?'] not in tfs[each[0]]:
                    tfs[each[0]].append([each[1], '?'])
            else:
                tfs[each[0]].append([each[1], direction_dict[each[2]]])
print("tf_dict Loaded")

for file in os.listdir(genie_dir):
    if '.csv' not in file:
        continue
    print('Working on: ', file)
    #genie = pd.read_csv('/Users/nathandmaulding/Desktop/EBI4_1_4.csv', index_col=0)
    genie = pd.read_csv(genie_dir + file, index_col=0)

    tfs2 = {}
    for k in tfs.keys():
        if k in genie['regulatoryGene'].tolist():
            tfs2[k] = tfs[k]

    # make random distributions
    d_lengths = [len(targets) for targets in tfs2.values()]
    mean_factors = (np.max(d_lengths))
    uTargets = genie.targetGene.unique()
    random.seed(11)
    rTargets1 = random.choices(uTargets, k=int(np.max(d_lengths)))
    random.seed(22)
    rTargets2 = random.choices(uTargets, k=int(np.max(d_lengths)))
    random.seed(33)
    rTargets3 = random.choices(uTargets, k=int(np.max(d_lengths)))
    background1 = [[x, '?'] for x in rTargets1]
    background2 = [[x, '?'] for x in rTargets2]
    background3 = [[x, '?'] for x in rTargets3]
    rbackground = [background1, background2, background3]
    random.seed(44)
    rFactors = random.choices(list(tfs2.keys()), k=3)
    randtfs = {rFactors[i]: rbackground[i] for i in range(len(rFactors))}
    rand_dist = []
    for tf in randtfs.keys():
        sub_genie = genie.loc[genie['regulatoryGene'] == tf]
        targets = []
        for t in randtfs[tf]:
            if t[0] not in targets:
                targets.append(t[0])
        tf_to_tar = sub_genie.loc[sub_genie['targetGene'].isin(targets)]
        weight_dist = tf_to_tar['weight'].tolist()
        rand_dist.extend(weight_dist)
    print('finished random')

    ks_dict = {}
    for tf in tfs2.keys():
        sub_genie = genie.loc[genie['regulatoryGene'] == tf]
        targets = []
        for t in tfs2[tf]:
            if t[0] not in targets:
                targets.append(t[0])
        tf_to_tar = sub_genie.loc[sub_genie['targetGene'].isin(targets)]
        weight_dist = tf_to_tar['weight'].tolist()
        if len(weight_dist) > 0:
            ks_dict[tf] = stats.ks_2samp(weight_dist, rand_dist)

    ks_df = pd.DataFrame.from_dict(ks_dict, orient='index')
    ks_df = ks_df.rename(index={0: "KSstat", 1: "pvalue"})
    ks_df = ks_df.sort_values('pvalue')

    suffix = '.csv'
    ks_df.to_csv(genie_dir + 'sig_factors/' + file[:file.rindex(suffix)] + '_genie_factors' + suffix)