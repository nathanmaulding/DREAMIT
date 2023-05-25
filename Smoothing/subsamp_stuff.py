import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

aicdf = pd.read_csv('/Users/nathandmaulding/Desktop/EBI5_2_3_aic_df.csv', index_col=0) #pd.read_csv('/Users/nathandmaulding/Desktop/ebi5_3_4_aic_df.csv', index_col=0)
print(aicdf)
aiccols = []
for col in aicdf.columns:
    aiccols.append(int(col))
aicdf.columns = aiccols

s = open('/Users/nathandmaulding/Desktop/EBI5_2_3_subsample_dictionary.txt', 'r') #open('/Users/nathandmaulding/Desktop/small_sub2.txt', 'r')
sub = ''
for i, line in enumerate(s):
    if i > 0:
        sub = sub + line
subdict = eval(sub)
print(subdict.keys())

def vanilla(datalofl):
    d = {}
    for elem in datalofl:
        for i in range(len(elem)):
            if i in d:
                d[i].append(elem[i])
            else:
                d[i] = [elem[i]]
    twosd = []
    ptwosd = []
    for k,v in d.items():
        #print(f'The average of {k} is {np.mean(v)} stdev of {np.std(v)} var of {np.std(v)*2} pvar {(np.std(v)*2)/abs(np.mean(v))}')
        if np.std(v) > 0:
            twosd.append(np.std(v))
            # coefficient of variation
            # https://en.wikipedia.org/wiki/Coefficient_of_variation
            ptwosd.append((np.std(v))/abs(np.mean(v)))
        else:
            twosd.append(0)
            ptwosd.append(0)
    return twosd, ptwosd

print(len(subdict.keys()))
g = list(subdict.keys())[0]
subdf = pd.DataFrame(subdict[g])
subcols = subdf.columns.tolist()
print(aiccols)
print(subcols)
assert aiccols == subcols

allcoeff = pd.DataFrame(columns=subdf.columns, index=subdf.index)
allknots = pd.DataFrame(columns=subdf.columns, index=subdf.index)
for gene in subdict.keys():
    subdf = pd.DataFrame(subdict[gene])
    dfcoeff = pd.DataFrame(columns=subdf.columns, index=subdf.index)
    dfknots = pd.DataFrame(columns=subdf.columns, index=subdf.index)
    for row in range(len(subdf.values)):
        for col in range(len(subdf.values[row])):
            sub_coeff = subdf.values[row][col][0]
            #print('subcoeff', len(sub_coeff), len(sub_coeff[0]), sub_coeff)
            coeff_twosd, pcoeff_twosd = vanilla(sub_coeff)
            #print('pcoeff_twosd', len(pcoeff_twosd), pcoeff_twosd)
            sub_knots = subdf.values[row][col][1]
            knot_twosd, pknot_twosd = vanilla(sub_knots)
            sub_res = subdf.values[row][col][2]
            res_var = np.std(sub_res)*3
            rowIndex = dfcoeff.index[row]
            colName = dfcoeff.columns[col]
            dfcoeff.loc[rowIndex, colName] = np.mean(pcoeff_twosd)
            dfknots.loc[rowIndex, colName] = np.mean(pknot_twosd)
    allcoeff = allcoeff.add(dfcoeff, fill_value=0)
    allknots = allknots.add(dfknots, fill_value=0)


allcoeff = allcoeff/len(subdict.keys())
allknots = allknots/len(subdict.keys())
print('allcoeff after', allcoeff)
print('allknots after', allknots)

opt_param_coeff = allcoeff.min().min()
opt_param_knots = allknots.min().min()
print('opt', opt_param_coeff, opt_param_knots)
for col in allcoeff.columns.tolist():
    if opt_param_coeff in allcoeff[col].tolist():
        row = allcoeff[col].tolist().index(opt_param_coeff)
        print('For coeff -- Smoothing Factor of ', allcoeff.index[row], 'Bin number of ', col)
    if opt_param_knots in allknots[col].tolist():
        row = allknots[col].tolist().index(opt_param_knots)
        print('For knots -- Smoothing Factor of ', allknots.index[row], 'Bin number of ', col)

print('coeff with aic -- ', aicdf.corrwith(allcoeff, axis=0))
print('coeff with aic -- ', aicdf.corrwith(allcoeff, axis=1))
print('\n')
print('coeff with aic -- ', aicdf.corrwith(allknots, axis=0))
print('coeff with aic -- ', aicdf.corrwith(allknots, axis=1))
print('\n')
print('coeff with knots -- ', allcoeff.corrwith(allknots, axis=0))
print('coeff with knots -- ', allcoeff.corrwith(allknots, axis=1))

allcoeff.to_csv('/Users/nathandmaulding/Desktop/EBI5_2_3_coeff_df.csv')
#allknots.to_csv('/Users/nathandmaulding/Desktop/ebi4_1_4_knots_df.csv')

distancedf = pd.DataFrame(columns=aicdf.columns, index=aicdf.index)
maxaic = aicdf.max().max()
maxcoeff = allcoeff.max().max()
for col in aicdf.columns:
    plt.scatter(aicdf[col].tolist(), allcoeff[col].tolist())
    for ind in aicdf.index:
        # normalize by the maximum point in the dataframe and calculate the distance from the origin
        # normalization allows for each axis to have an equal effect and distance is used to balance the prefential treatment of one component over another
        distancedf[col][ind] = np.sqrt(((aicdf[col][ind]/maxaic)**2) + ((allcoeff[col][ind]/maxcoeff)**2)) #np.sqrt( ((x2 - 0)**2) + ((y2 - 0)**2) )
print(distancedf)
best = distancedf.min().min()
worst = distancedf.max().max()
smoothing_factors = [0.01, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 7.5, 10]
for col in distancedf.columns.tolist():
    if best in distancedf[col].tolist():
        row = distancedf[col].tolist().index(best)
        print('best Smoothing Factor of ', smoothing_factors[row], 'Bin number of ', col)
    if worst in distancedf[col].tolist():
        row = distancedf[col].tolist().index(worst)
        print('worst Smoothing Factor of ', smoothing_factors[row], 'Bin number of ', col)


#plt.ylim([0,50])
plt.xlabel('AIC')
plt.ylabel('Coefficient of Variation')
plt.savefig('AICvsCoeff_EBI4_4_3.png', dpi=300)
plt.show()
