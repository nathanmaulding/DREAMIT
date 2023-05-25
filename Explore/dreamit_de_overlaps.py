import sys
import os
import pandas as pd


#folder = '/Users/nathandmaulding/Desktop/DREAMIT2/Data/Friedman_v3_bins20/'
#folder_de = folder + 'diffexp/'
#folder_dreamit = folder + 'new_dreamit_results/'
#folder_path = folder + 'pathway/'
folder = '/Users/nathandmaulding/Desktop/v20Data/'
sub = 'friedman/'
folder_de = folder + sub
folder_dreamit = folder + sub
#folder_path = folder + 'pathway/'
tissue = 'Heart'
cell_type = 'Stem cell'

factor_markers = []
for line in open('/Users/nathandmaulding/Desktop/DREAMIT/TF-Marker-All-TF-markers.txt', 'r'):
    new = line.split('\t')
    if len(new) > 6:
        if new[5] == tissue:
        #if new[3] == cell_type:
            if new[1] not in factor_markers:
                factor_markers.append(new[1])
#factor_markers = [i.title() for i in factor_markers]  #### this is only for PAUL data
# also check lines 47-50
print('TF Markers of ' + tissue + ':', len(factor_markers), factor_markers)

def openFile(path):
    data = open(path, 'r')
    figDict = ''
    branch = ''
    for i, line in enumerate(data):
        if i == 0:
            branch = line.strip('\n')
        if i >= 3:
            figDict = figDict + line
    figDict = eval(figDict)
    return figDict, branch

fs = []
tars = []
for each in open('/Users/nathandmaulding/Desktop/DREAMIT/trrust_rawdata.human.tsv', 'r'):
    e = each.split()
    ### .title() for PAUL data
    if e[0].title() not in fs:
        fs.append(e[0].title())
    if e[1].title() not in tars:
        tars.append(e[1].title())


for file in os.listdir(folder_de):
    if 'Smoothed_Expression_dictionary' in file:
        path = folder_de + file
        exp, branch = openFile(path)
        ftested = []
        ttested = []
        for each_gene in exp.keys():
            if sum(exp[each_gene]) > 0:
                if each_gene in fs and each_gene not in ftested:
                    ftested.append(each_gene)
                if each_gene in tars and each_gene not in ttested:
                    ttested.append(each_gene)
        marktested = [i for i in ftested if i in factor_markers]
        print("Tested on Branch: ", branch, "Factors tested: ", len(ftested), "Targets tested: ", len(ttested), "Markers tested: ", len(marktested))

def diffExpReport(filename):
    de = {}
    for file in os.listdir(folder_de):
        if filename in file:
            df = pd.read_csv(folder_de + file, header=2, index_col=0, sep='\t')
            sig = df[df['FDR'] <= 0.05]
            fsig = sig[sig.index.isin(ftested)]
            factors = fsig.shape[0]
            tsig = sig[sig.index.isin(ttested)]
            targets = tsig.shape[0]
            branch = file.split('_' + filename + '.txt')[0]
            de[branch] = fsig
            de_marks = [i for i in fsig.index if i in factor_markers]
            print(filename, " on Branch: ", branch, "Sig Factors: ", factors, "Sig Targets: ", targets, "Sig Markers: ", len(de_marks))
    return de

raw_de = diffExpReport('raw_diffexp')
smooth_de = diffExpReport('smooth_diffexp')

dreamit = {}
for file in os.listdir(folder_dreamit):
    if 'Summary_of_Significantly_NonRandom_Findings' in file:
        df2 = pd.read_csv(folder_dreamit + file, header=2, index_col=0, sep='\t')
        if len(df2["Branch"]) > 0:
            branch = file.split('_Summary_of_Significantly_NonRandom_Findings.txt')[0]
            dreamit[branch] = df2
            tars = []
            for i in range(len(df2["Names of Targets"])):
                targets = df2["Names of Targets"][i]
                for each in targets:
                    if each not in tars:
                        tars.append(each)
            ks05 = df2[df2['KS Statistic Value'] >= 0.5]
            ks075 = df2[df2['KS Statistic Value'] >= 0.75]
            marks = [i for i in df2.index.unique() if i in factor_markers]
            marks05 = [i for i in ks05.index.unique() if i in factor_markers]
            marks075 = [i for i in ks075.index.unique() if i in factor_markers]
            print("DREAMIT2 on branch: ", branch, "Total Sig Factors: ", len(df2.index), "Unique Sig Factors: ", len(df2.index.unique()), "Unique Targets: ", len(tars))
            print("DREAMIT2 on branch: ", branch, "KS0 Factors: ", len(df2.index.unique()), "KS0.5 Factors: ", len(ks05.index.unique()), "KS0.75 Factors: ", len(ks075.index.unique()))
            print("DREAMIT2 on branch: ", branch, "Marker KS0: ", len(marks), "Marker KS0.5: ", len(marks05), "Marker KS0.75: ", len(marks075))
#print("de", de)
only = {}
overlap_smooth = {}
overlap_raw = {}
for branch in dreamit.keys():
    smoothde_factors = smooth_de[branch].index.tolist()
    rawde_factors = raw_de[branch].index.tolist()
    dreamit_factors = dreamit[branch].index.unique().tolist()

    smooth_both = [x for x in smoothde_factors if x in dreamit_factors]
    onlysmoothde = [x for x in smoothde_factors if x not in dreamit_factors]
    nosmooth_dreamit = [x for x in dreamit_factors if x not in smoothde_factors]
    raw_both = [x for x in rawde_factors if x in dreamit_factors]
    onlyrawde = [x for x in rawde_factors if x not in dreamit_factors]
    noraw_dreamit = [x for x in dreamit_factors if x not in rawde_factors]
    no_raworsmooth_dreamit = [x for x in dreamit_factors if x not in rawde_factors and x not in smoothde_factors]
    smooth_both_markers = [x for x in smoothde_factors if x in dreamit_factors and x in factor_markers]
    onlysmoothde_markers = [x for x in smoothde_factors if x not in dreamit_factors and x in factor_markers]
    nosmooth_dreamit_markers = [x for x in dreamit_factors if x not in smoothde_factors and x in factor_markers]
    raw_both_markers = [x for x in rawde_factors if x in dreamit_factors and x in factor_markers]
    onlyrawde_markers = [x for x in rawde_factors if x not in dreamit_factors and x in factor_markers]
    noraw_dreamit_markers = [x for x in dreamit_factors if x not in rawde_factors and x in factor_markers]
    no_raworsmooth_dreamit_markers = [x for x in dreamit_factors if x not in rawde_factors and x not in smoothde_factors and x in factor_markers]
    print("Branch: ", branch, "Factors found in smoothDE and DREAMIT2: ", len(smooth_both), smooth_both)
    print("Branch: ", branch, "Factors only in smoothDE: ", len(onlysmoothde), onlysmoothde)
    print("Branch: ", branch, "Factors found in rawDE and DREAMIT2: ", len(raw_both), raw_both)
    print("Branch: ", branch, "Factors only in rawDE: ", len(onlyrawde), onlyrawde)
    print("Branch: ", branch, "Factors missed by smoothDE and found by DREAMIT2: ", len(nosmooth_dreamit), nosmooth_dreamit)
    print("Branch: ", branch, "Factors missed by rawDE and found by DREAMIT2: ", len(noraw_dreamit), noraw_dreamit)
    print("Branch: ", branch, "Factors missed by smoothDE and rawDE and found by DREAMIT2: ", len(no_raworsmooth_dreamit), no_raworsmooth_dreamit)
    print("Branch: ", branch, "Markers found in smoothDE and DREAMIT2: ", len(smooth_both_markers), smooth_both_markers)
    print("Branch: ", branch, "Markers only in smoothDE: ", len(onlysmoothde_markers), onlysmoothde_markers)
    print("Branch: ", branch, "Markers found in rawDE and DREAMIT2: ", len(raw_both_markers), raw_both_markers)
    print("Branch: ", branch, "Markers only in rawDE: ", len(onlyrawde_markers), onlyrawde_markers)
    print("Branch: ", branch, "Markers missed by smoothDE and found by DREAMIT2: ", len(nosmooth_dreamit_markers), nosmooth_dreamit_markers)
    print("Branch: ", branch, "Markers missed by rawDE and found by DREAMIT2: ", len(noraw_dreamit_markers), noraw_dreamit_markers)
    print("Branch: ", branch, "Markers missed by smoothDE and rawDE and found by DREAMIT2: ", len(no_raworsmooth_dreamit_markers), no_raworsmooth_dreamit_markers)
    print('\n')
    df = dreamit[branch]
    df_onlydreamit = df[df.index.isin(no_raworsmooth_dreamit)]
    df_smoothboth = df[df.index.isin(smooth_both)]
    df_rawboth = df[df.index.isin(raw_both)]
    only[branch] = {}
    overlap_smooth[branch] = {}
    overlap_raw[branch] = {}
    for each in nosmooth_dreamit:
        only[branch][each] = df_onlydreamit[df_onlydreamit.index == each]["Method"].to_list()
    for factor in smooth_both:
        overlap_smooth[branch][factor] = df_smoothboth[df_smoothboth.index == factor]["Method"].to_list()
    for f in raw_both:
        overlap_raw[branch][f] = df_rawboth[df_rawboth.index == f]["Method"].tolist()

for branch in only.keys():
    print('\n')
    print("ONLY DREAMIT2 FACTORS AND METHODS: ", branch)
    for f in only[branch].keys():
        print("Factor: ", f, "Methods: ", only[branch][f])

# Overlap between each method and DE
pear, allpear = 0, 0
spear, allspear = 0, 0
dtw, alldtw = 0, 0
mi, allmi = 0, 0
roll, allroll = 0, 0
for branch in overlap_smooth.keys():
    all = dreamit[branch]['Method'].tolist()
    de_factors = smooth_de[branch].index.tolist()
    for methods in overlap_smooth[branch].values():
        for m in methods:
            if 'Pearson' in m and 'Rolling' not in m:
                pear += 1
            if 'Spearman' in m and 'Rolling' not in m:
                spear += 1
            if 'Dynamic' in m and 'Rolling' not in m:
                dtw += 1
            if 'Mutual' in m and 'Rolling' not in m:
                mi += 1
            if 'Rolling' in m:
                roll += 1
    for each in all:
        if 'Pearson' in each and 'Rolling' not in each:
            allpear += 1
        if 'Spearman' in each and 'Rolling' not in each:
            allspear += 1
        if 'Dynamic' in each and 'Rolling' not in each:
            alldtw += 1
        if 'Mutual' in each and 'Rolling' not in each:
            allmi += 1
        if 'Rolling' in each:
            allroll += 1
    print('\n')
    print("OVERLAP BETWEEN DREAMIT2 PEARSON AND smoothDE:", branch, " -- # of overlapping factors:", pear, '-- # of all found by method:', allpear, "-- # of all de:", len(de_factors))
    print("OVERLAP BETWEEN DREAMIT2 SPEARMAN AND smoothDE:", branch, " -- # of overlapping factors:", spear, '-- # of all found by method:', allspear, "-- # of all de:", len(de_factors))
    print("OVERLAP BETWEEN DREAMIT2 DTW AND smoothDE:", branch, " -- # of overlapping factors:", dtw, '-- # of all found by method:', alldtw, "-- # of all de:", len(de_factors))
    print("OVERLAP BETWEEN DREAMIT2 MUTUAL INFORMATION AND smoothDE:", branch, " -- # of overlapping factors:", mi, '-- # of all found by method:', allmi, "-- # of all de:", len(de_factors))
    print("OVERLAP BETWEEN DREAMIT2 ROLLING AND smoothDE:", branch, " -- # of overlapping factors:", roll, '-- # of all found by method:', allroll, "-- # of all de:", len(de_factors))

# Overlap between each method and rawDE
pear, allpear = 0, 0
spear, allspear = 0, 0
dtw, alldtw = 0, 0
mi, allmi = 0, 0
roll, allroll = 0, 0
for branch in overlap_raw.keys():
    all = dreamit[branch]['Method'].tolist()
    de_factors = raw_de[branch].index.tolist()
    for methods in overlap_raw[branch].values():
        for m in methods:
            if 'Pearson' in m and 'Rolling' not in m:
                pear += 1
            if 'Spearman' in m and 'Rolling' not in m:
                spear += 1
            if 'Dynamic' in m and 'Rolling' not in m:
                dtw += 1
            if 'Mutual' in m and 'Rolling' not in m:
                mi += 1
            if 'Rolling' in m:
                roll += 1
    for each in all:
        if 'Pearson' in each and 'Rolling' not in each:
            allpear += 1
        if 'Spearman' in each and 'Rolling' not in each:
            allspear += 1
        if 'Dynamic' in each and 'Rolling' not in each:
            alldtw += 1
        if 'Mutual' in each and 'Rolling' not in each:
            allmi += 1
        if 'Rolling' in each:
            allroll += 1
    print('\n')
    print("OVERLAP BETWEEN DREAMIT2 PEARSON AND rawDE:", branch, " -- # of overlapping factors:", pear, '-- # of all found by method:', allpear, "-- # of all de:", len(de_factors))
    print("OVERLAP BETWEEN DREAMIT2 SPEARMAN AND rawDE:", branch, " -- # of overlapping factors:", spear, '-- # of all found by method:', allspear, "-- # of all de:", len(de_factors))
    print("OVERLAP BETWEEN DREAMIT2 DTW AND rawDE:", branch, " -- # of overlapping factors:", dtw, '-- # of all found by method:', alldtw, "-- # of all de:", len(de_factors))
    print("OVERLAP BETWEEN DREAMIT2 MUTUAL INFORMATION AND rawDE:", branch, " -- # of overlapping factors:", mi, '-- # of all found by method:', allmi, "-- # of all de:", len(de_factors))
    print("OVERLAP BETWEEN DREAMIT2 ROLLING AND rawDE:", branch, " -- # of overlapping factors:", roll, '-- # of all found by method:', allroll, "-- # of all de:", len(de_factors))

# the following is for pathways
'''
direct = ['pearson', 'spearman']
for file in os.listdir(folder_path):
    if 'modules' in file:
        df = pd.read_csv(folder_path + file, header=2, index_col=0, sep='\t')
        branch = file.split('_modules.txt')[0]
        ns = [i for i in range(len(df.index.tolist())) if df.index.tolist()[i] == 'Not Significant']
        df = df[:ns[0]]
        df['FDR'] = pd.to_numeric(df['FDR'])
        df005 = df[df['FDR'] <= 0.05]
        u005 = df005.index.unique().tolist()
        df001 = df[df['FDR'] <= 0.01]
        u001 = df001.index.unique().tolist()
        df_direct = df005[df005['Method'].isin(direct)]
        udirect = df_direct.index.unique().tolist()
        print("Pathways on Branch: ", branch, "UPaths FDR005: ", len(u005), "UPaths FDR001: ", len(u001), "UPaths Direct: ", len(udirect))
'''
