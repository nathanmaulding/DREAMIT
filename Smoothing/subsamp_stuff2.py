# This will plot the distance of each param selection in terms of CV and AIC

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

aicdf = pd.read_csv('/Users/nathandmaulding/Desktop/splinechoice/paul_withregs_stem_to_erythro_aic_df.csv', index_col=0)
print(aicdf)
s = open('/Users/nathandmaulding/Desktop/splinechoice/paul_withregs_stem_to_erythro_subsample_dictionary.txt', 'r') #open('/Users/nathandmaulding/Desktop/small_sub2.txt', 'r')
sub = ''
for i, line in enumerate(s):
    if i > 0:
        sub = sub + line
subdict = eval(sub)
print(subdict.keys())

allcoeff = pd.DataFrame(columns=aicdf.columns, index=aicdf.index)
#allcoeff = pd.DataFrame()
for gene, bin_dicts in subdict.items():
    dfcoeff = pd.DataFrame(columns=aicdf.columns, index=aicdf.index)
    for bin, sf_dicts in bin_dicts.items():
        sf_coeff = []
        for sf, spl_attr in sf_dicts.items():
            d = {}
            pcoeffvar = []
            for elem in spl_attr[0]:
                for i in range(len(elem)):
                    if i in d:
                        d[i].append(elem[i])
                    else:
                        d[i] = [elem[i]]
            for k, v in d.items():
                if np.std(v) > 0:
                    pcoeffvar.append((np.std(v)) / abs(np.mean(v)))
                else:
                    pcoeffvar.append(0)
            sf_coeff.append(np.mean(pcoeffvar))
        dfcoeff[str(bin)] = sf_coeff
    allcoeff = allcoeff.add(dfcoeff, fill_value=0)

#print(subdict.keys())
print('allcoeff')
#allcoeff.to_csv('/Users/nathandmaulding/Desktop/paul_withregs_stem_to_mono_coeff_df.csv')
print(allcoeff.columns.tolist())
#print(allcoeff['4'], allcoeff[4])
print(aicdf.columns.tolist())
#allcoeff = pd.read_csv('/Users/nathandmaulding/Desktop/subsample/data/ebi5_2_3_coeff_df.csv', index_col=0)
#allknots = pd.read_csv('/Users/nathandmaulding/Desktop/subsample/data/ebi5_2_3_knots_df.csv', index_col=0)



'''
print('coeff with aic -- ', aicdf.corrwith(allcoeff, axis=0))
print('coeff with aic -- ', aicdf.corrwith(allcoeff, axis=1))
print('\n')
print('knots with aic -- ', aicdf.corrwith(allknots, axis=0))
print('knots with aic -- ', aicdf.corrwith(allknots, axis=1))
print('\n')
print('coeff with knots -- ', allcoeff.corrwith(allknots, axis=0))
print('coeff with knots -- ', allcoeff.corrwith(allknots, axis=1))
'''


aicdf = aicdf.div(100)
allcoeff = allcoeff.div(100)
distancedf = pd.DataFrame(columns=aicdf.columns, index=aicdf.index)
distancedf2 = pd.DataFrame(columns=aicdf.columns, index=aicdf.index)
maxaic = aicdf.max().max()
maxcv = allcoeff.max().max()
coeff_thresh = 1
for col in aicdf.columns:
    for ind in aicdf.index:
        # normalize by the maximum point in the dataframe and calculate the distance from the origin
        # normalization allows for each axis to have an equal effect and distance is used to balance the prefential treatment of one component over another
        distancedf[col][ind] = np.sqrt(((aicdf[col][ind] / maxaic) ** 2) + ((allcoeff[col][ind] / maxcv) ** 2))  # np.sqrt( ((x2 - 0)**2) + ((y2 - 0)**2) )
        if allcoeff[col][ind] <= coeff_thresh:
            distancedf2[col][ind] = np.sqrt(((aicdf[col][ind] / maxaic) ** 2) + ((allcoeff[col][ind] / maxcv) ** 2))  # np.sqrt( ((x2 - 0)**2) + ((y2 - 0)**2) )
print('d', distancedf)
print('d2', distancedf2)

best_param = distancedf.min().min()
worst_param = distancedf.max().max()
for col in distancedf.columns.tolist():
    if best_param in distancedf[col].tolist() and best_param >= 0:
        row = distancedf[col].tolist().index(best_param)
        print('Best - ', 'Smoothing Factor of ', distancedf.index[row], 'Bin number of ', col)
    if worst_param in distancedf[col].tolist() and worst_param >= 0:
        row = distancedf[col].tolist().index(worst_param)
        print('Worst - ', 'Smoothing Factor of ', distancedf.index[row], 'Bin number of ', col)

best_param_cv = allcoeff.min().min()
worst_param_cv = allcoeff.max().max()
for col in allcoeff.columns.tolist():
    if best_param_cv in allcoeff[col].tolist() and best_param_cv >= 0:
        row = allcoeff[col].tolist().index(best_param_cv)
        print('Best Coeff- ', 'Smoothing Factor of ', allcoeff.index[row], 'Bin number of ', col)
    if worst_param_cv in allcoeff[col].tolist() and worst_param_cv >= 0:
        row = allcoeff[col].tolist().index(worst_param_cv)
        print('Worst Coeff- ', 'Smoothing Factor of ', allcoeff.index[row], 'Bin number of ', col)


labels = [i + ' bins' for i in aicdf.columns.tolist()]
labels += [str(i) + ' SF' for i in aicdf.index.tolist()]
colors = ['black', 'blue', 'red', 'green', 'yellow', 'pink', 'darkviolet', 'teal']
markers = ["*", "^", "o", 'v', 's', 'p', 'd', 'x', '+', 'P']
for i, col in enumerate(aicdf.columns):
    for j, ind in enumerate(aicdf.index):
        plt.scatter(aicdf[col][ind], allcoeff[col][ind], color=colors[i], marker=markers[j])

f = lambda m,c: plt.plot([],[],marker=m, color=c, ls="none")[0]
handles = [f("s", colors[i]) for i in range(len(colors))]
handles += [f(markers[i], "k") for i in range(len(markers))]

plt.legend(handles, labels, loc='upper right', fontsize=6)
plt.ylim([0, 1])
#plt.xlim([0, (maxaic/100)+1])
plt.xlabel('Akaike Information Criterion (AIC)')
plt.ylabel('Coefficient of Variation (CV)')
plt.title('AIC vs CV with threshold at CV > 1')
plt.savefig('AICvsCoeff_Paulwithregs_erythro_thresh.png', dpi=300)
plt.show()


'''
# Creating figure
fig = plt.figure(figsize=(10, 7))
ax = plt.axes(projection="3d")

# Creating plot
ax.scatter3D(aicdf, allcoeff, allknots, color="green")
plt.title("simple 3D scatter plot")
ax.set_ylim([0,20])
ax.set_xlabel('AIC', fontweight ='bold')
ax.set_ylabel('Coeff', fontweight ='bold')
ax.set_zlabel('Knots', fontweight ='bold')

# show plot
plt.show()

### FIND POINT THAT IS CLOSEST TO ORIGIN (AIC=0, Coeff=0) or (AIC=0, Coeff=0, Knots=0) ###
'''
