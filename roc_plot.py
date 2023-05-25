import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import os

# FDR (x-axis) vs Positive Rate (y-axis)
# should be done for all tested factors
# and for all factor markers
folder = '/Users/nathandmaulding/Desktop/DREAMIT2/Data/'

'''
# for roc_data.py
ldreamit = []
for each in open('/Users/nathandmaulding/Desktop/dreamit_fdr_results.txt', 'r'):
    dreamit = each.split()
    for i in dreamit:
        x = float(i)
        ldreamit.append(x)
print(ldreamit)

lde = []
for each in open('/Users/nathandmaulding/Desktop/de_fdr_results.txt', 'r'):
    de = each.split()
    for i in de:
        x = float(i)
        lde.append(x)
print(lde)
'''

# for roc_data2.py
dreamit_df = pd.read_csv('/Users/nathandmaulding/Desktop/dreamit_fdr_results2.csv', index_col=0)
de_df = pd.read_csv('/Users/nathandmaulding/Desktop/de_fdr_results2.csv', index_col=0)


#tissue = 'Brain'
#species = 'Human'
#species_dict = {'EBI': 'Human', 'friedman': 'Human', 'paul': 'Mouse'}
tissue_dict = {'EBI3': 'Testis', 'EBI4': 'Brain', 'EBI5': 'Bone marrow', 'EBI6': 'Fetal liver', 'EBI7': 'Embryo', 'EBI8': 'Embryo', 'EBI9': 'Retina', 'EBI10': 'Brain', 'friedman': 'Heart', 'paul': 'Embryo'}

'''
# make dataset
lde = []
ldreamit = []
for each in tissue_dict.keys():
    for col in dreamit_df.columns.tolist():
        if each in col:
            if 'EBI' in col or 'friedman' in col:
                species = 'Human'
            else:
                species = 'Mouse'
            tissue = tissue_dict[each]
            markers = {}
            for line in open('/Users/nathandmaulding/Desktop/DREAMIT2/Single_cell_markers.txt', 'r'):
                new = line.split('\t')
                if new[1] == tissue and new[0] == species:
                    genes = new[8].split(', ')
                    for i in genes:
                        if i not in markers:
                            markers[i] = [new[5]]
                        else:
                            markers[i].append(new[5])
            factor_markers = markers.keys()
            sub_de_df = de_df[de_df.index.isin(factor_markers)]
            if sub_de_df.empty:
                print('de', tissue, factor_markers)
            sub_dreamit_df = dreamit_df[dreamit_df.index.isin(factor_markers)]
            if sub_dreamit_df.empty:
                print('dreamit', tissue, factor_markers)
            ldreamit = ldreamit + sub_dreamit_df[col].dropna().tolist()
            lde = lde + sub_de_df[col].dropna().tolist()
print(len(ldreamit), ldreamit)
print(len(lde), lde)
'''

'''
# make dataset
lde = []
lde2 = []
ldreamit = []
ldreamit_full = []
for each in tissue_dict.keys():
    for col in dreamit_df.columns.tolist():
        if each in col:
            if 'EBI' in col or 'friedman' in col:
                species = 'Human'
            else:
                species = 'Mouse'
            tissue = tissue_dict[each]
            markers = {}
            for line in open('/Users/nathandmaulding/Desktop/DREAMIT2/Single_cell_markers.txt', 'r'):
                new = line.split('\t')
                if new[1] == tissue and new[0] == species:
                    genes = new[8].split(', ')
                    for i in genes:
                        if i not in markers:
                            markers[i] = [new[5]]
                        else:
                            markers[i].append(new[5])
            factor_markers = markers.keys()
            sub_de_df = de_df[de_df.index.isin(factor_markers)]
            if len(sub_de_df.index) > 0:
                #print('de', len(de_df.index), len(sub_de_df.index))
                lde2 = lde2 + de_df[col].dropna().tolist()
            sub_dreamit_df = dreamit_df[dreamit_df.index.isin(factor_markers)]
            if len(sub_dreamit_df.index) > 0:
                #print('dreamit', len(dreamit_df.index), len(sub_dreamit_df.index))
                ldreamit_full = ldreamit_full + dreamit_df[col].dropna().tolist()
            ldreamit = ldreamit + sub_dreamit_df[col].dropna().tolist()
            lde = lde + sub_de_df[col].dropna().tolist()
print(len(ldreamit), ldreamit)
print(len(ldreamit_full), ldreamit_full)
print(len(lde), lde)
print(len(lde2), lde2)
'''

# make dataset
lde = []
lde2 = []
ldreamit = []
ldreamit_full = []
for each in tissue_dict.keys():
    for col in dreamit_df.columns.tolist():
        if each in col:
            if 'EBI' in col or 'friedman' in col:
                species = 'Human'
            else:
                species = 'Mouse'
            tissue = tissue_dict[each]
            markers = {}
            for line in open('/Users/nathandmaulding/Desktop/DREAMIT/Single_cell_markers.txt', 'r'):
                new = line.split('\t')
                if new[1] == tissue and new[0] == species:
                    genes = new[8].split(', ')
                    for i in genes:
                        if i not in markers:
                            markers[i] = [new[5]]
                        else:
                            markers[i].append(new[5])
            factor_markers = markers.keys()
            sub_de_df = de_df[de_df.index.isin(factor_markers)]
            if len(sub_de_df[col].dropna().index) > 0:
                #print('de', len(de_df.index), len(sub_de_df.index))
                print(tissue, factor_markers)
                print('de', de_df[col].dropna())
                print(tissue, len(factor_markers), factor_markers)
                print('sub_de', sub_de_df[col].dropna())
                lde2.append(de_df[col].dropna().tolist())
                lde.append(sub_de_df[col].dropna().tolist())

            sub_dreamit_df = dreamit_df[dreamit_df.index.isin(factor_markers)]
            if len(sub_dreamit_df[col].dropna().index) > 0:
                #print('dreamit', len(dreamit_df.index), len(sub_dreamit_df.index))
                print(tissue, factor_markers)
                print('dreamit', dreamit_df[col].dropna())
                print(tissue, len(factor_markers), factor_markers)
                print('sub_dreamit', sub_dreamit_df[col].dropna())
                ldreamit_full.append(dreamit_df[col].dropna().tolist())
                ldreamit.append(sub_dreamit_df[col].dropna().tolist())

print(len(ldreamit), ldreamit)  # this is the marker genes
print(len(ldreamit_full), ldreamit_full)  # this is the total genes
print(len(lde), lde)  # this is the marker genes
print(len(lde2), lde2)  # this is the total genes

# figure specs
figureHeight = 6
figureWidth = 12
plt.figure(figsize=(figureWidth, figureHeight))
panelWidth = 8
panelHeight = 4
relPanWidth = panelWidth/figureWidth
relPanHeight = panelHeight/figureHeight

# main panel specs
panel1 = plt.axes([0.09, 0.25, relPanWidth, relPanHeight])  # 0.04, 0.15
panel1.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True, right=False, labelright=False, top=False, labeltop=False)
#panel2 = plt.axes([0.09, 0.25, relPanWidth, relPanHeight])  # 0.04, 0.15
#panel2.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True, right=False, labelright=False, top=False, labeltop=False)



# plotting
#vals = [0, 1e-50, 1e-40, 1e-30, 1e-20, 1e-15, 1e-10, 1e-5, 1e-4, 1e-3, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1]
vals = [0, 1e-50, 1e-40, 1e-30, 1e-20, 1e-15, 1e-10, 1e-5, 1e-4, 1e-3, 0.01, 0.025, 0.05]
dreamit_ys = []
dreamit_ys2 = []
de_ys = []
de_ys2 = []
for j in range(len(ldreamit)):
    print('dreamit markers', j, len([x for x in ldreamit[j] if x <= 0.05])/len(ldreamit[j]), ldreamit[j])
    if len([x for x in ldreamit_full[j] if x <= 0.05]) > 0:
        print('dreamit total', len([x for x in ldreamit[j] if x <= 0.05]) / len([x for x in ldreamit_full[j] if x <= 0.05]))
    else:
        print('dreamit total', 0)
    print('de markers', len([x for x in lde[j] if x <= 0.05])/len(lde[j]))
    if len([x for x in lde2[j] if x <= 0.05]) > 0:
        print('de total', len([x for x in lde[j] if x <= 0.05]) / len([x for x in lde2[j] if x <= 0.05]))
    else:
        print('de total', 0)
    print()
    dreamit_ys.append([])
    dreamit_ys2.append([])
    de_ys.append([])
    de_ys2.append([])
    for i in vals:
        dreamit_ys[j].append(len([x for x in ldreamit[j] if x <= i])/len(ldreamit[j]))
        if len([x for x in ldreamit_full[j] if x <= i]) > 0:
            dreamit_ys2[j].append(len([x for x in ldreamit[j] if x <= i]) / len([x for x in ldreamit_full[j] if x <= i]))
        else:
            dreamit_ys2[j].append(0)
        de_ys[j].append(len([x for x in lde[j] if x <= i])/len(lde[j]))
        if len([x for x in lde2[j] if x <= i]) > 0:
            de_ys2[j].append(len([x for x in lde[j] if x <= i]) / len([x for x in lde2[j] if x <= i]))
        else:
            de_ys2[j].append(0)
print(dreamit_ys2, de_ys2)

average_dreamit = []
average_de = []
for i in range(len(dreamit_ys2[0])):
    dreamit_x = np.mean([j[i] for j in dreamit_ys2])
    de_x = np.mean([j[i] for j in de_ys2])
    average_dreamit.append(dreamit_x)
    average_de.append(de_x)

#################################
# all only
for x in range(len(dreamit_ys2)):

    #plt.scatter(vals, dreamit_ys2[x], alpha=0.5, color='black')
    #plt.scatter(vals, de_ys2[x], alpha=0.5, color='red')
    plt.plot(vals, dreamit_ys2[x], linestyle='solid', label='DREAMIT_' + str(x), alpha=0.5, color='black')
    plt.plot(vals, de_ys2[x], linestyle='dashed', label='DE_' + str(x), alpha=0.5, color='red')

    plt.legend(loc='lower right', fontsize=8)
    panel1.set_title('Markers found out of Total Results ' + str(x))
    panel1.set_xlabel('FDR')
    panel1.set_ylabel('Fraction of Markers in Results ' + str(x))

    plt.savefig('roc4_' + str(x) + '.png', dpi=300)
    plt.show()


###############################
# average and all
for x in range(len(dreamit_ys2)):

    plt.plot(vals, dreamit_ys2[x], linestyle='solid', alpha=0.1, color='black')
    plt.plot(vals, de_ys2[x], linestyle='dashed', alpha=0.1, color='red')

plt.plot(vals, average_dreamit, linestyle='solid', label='DREAMIT_mean', alpha=1, color='black')
plt.plot(vals, average_de, linestyle='dashed', label='DE_mean', alpha=1, color='red')
plt.legend(loc='lower right', fontsize=8)
plt.title('Markers found out of Total Results')
plt.xlabel('FDR')
plt.ylabel('Fraction of Markers in Results')

plt.savefig('roc4_average_and_all.png', dpi=300)
plt.show()

#################
# only average
plt.plot(vals, average_dreamit, linestyle='solid', label='DREAMIT_mean', alpha=1, color='black')
plt.plot(vals, average_de, linestyle='dashed', label='DE_mean', alpha=1, color='red')
plt.legend(loc='lower right', fontsize=8)
plt.title('Markers found out of Total Results')
plt.xlabel('FDR')
plt.ylabel('Fraction of Markers in Results')

plt.savefig('roc4_average_only.png', dpi=300)
plt.show()