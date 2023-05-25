import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import os

# for roc_data2.py
dreamit_df = pd.read_csv('/Users/nathandmaulding/Desktop/specificity_results/v20/fdr/dreamit_fdr_results_v20.csv',index_col=0)  # old data: dreamit_fdr_results3.csv
de_df = pd.read_csv('/Users/nathandmaulding/Desktop/specificity_results/v20/fdr/TargetRawDE_fdr_results.csv',index_col=0)  # old data: de_fdr_results3.csv
genie_df = pd.read_csv('/Users/nathandmaulding/Desktop/specificity_results/v20/fdr/TargetSmoothDE_fdr_results.csv',index_col=0)
output = '/Users/nathandmaulding/Desktop/specificity_results/v20/scratch/'
label1 = 'DREAMIT'
label2 = 'rawDEtargets'
label3 = 'smoothDEtargets'

#dreamit_df.columns = [i.replace('_v20', '') for i in dreamit_df.columns.tolist()]
#de_df.columns = [i.replace('_genie_pathways.csv', '') for i in de_df.columns.tolist()]
#genie_df.columns = [i.replace('_genie_pathways.csv', '') for i in genie_df.columns.tolist()]

tissue_dict = {'EBI3': 'Testis', 'EBI4': 'Brain', 'EBI5': 'Bone marrow', 'EBI6': 'Liver', 'EBI7': 'Embryo',
               'EBI8': 'Embryo', 'EBI9': 'Retina', 'EBI10': 'Brain', 'friedman': 'Heart', 'paul': 'Blood'}

de_total = []
de_total_positive = []
de_total_negative = []
dreamit_total = []
dreamit_total_positive = []
dreamit_total_negative = []
genie_total = []
genie_total_positive = []
genie_total_negative = []
data_details = []
de_total_factors = []
dreamit_total_factors = []
genie_total_factors = []
for each in tissue_dict.keys():
    for col in dreamit_df.columns.tolist():
        if each in col:
            coldreamit = dreamit_df[col].dropna()
            colde = de_df[col].dropna()
            colgenie = genie_df[col].dropna()
            if 'EBI' in col or 'friedman' in col:
                species = 'Human'
            else:
                species = 'Mouse'
            tissue = tissue_dict[each]
            factor_markers = []
            for line in open('/Users/nathandmaulding/Desktop/DREAMIT/TF-Marker-All-TF-markers.txt', 'r'):
                new = line.split('\t')
                if len(new) > 6:
                    if new[5] == tissue:
                        if new[1] not in factor_markers:
                            factor_markers.append(new[1])
            print('Tissue: ', tissue, factor_markers)
            markers2 = []
            omit_markers = []
            for fact in factor_markers:
                if fact in coldreamit.index.tolist() and fact in colde.index.tolist() and fact in colgenie.index.tolist():
                    markers2.append(fact)
                else:
                    omit_markers.append(fact)
            dreamit_negative = [x for x in coldreamit.index.tolist() if x not in factor_markers]
            de_negative = [x for x in colde.index.tolist() if x not in factor_markers]
            genie_negative = [x for x in colgenie.index.tolist() if x not in factor_markers]
            all_negative = set(dreamit_negative).intersection(de_negative, genie_negative)

            coldreamit_positive = coldreamit[coldreamit.index.isin(markers2)]
            coldreamit_negative = coldreamit[coldreamit.index.isin(all_negative)]
            colde_positive = colde[colde.index.isin(markers2)]
            colde_negative = colde[colde.index.isin(all_negative)]
            colgenie_positive = colgenie[colgenie.index.isin(markers2)]
            colgenie_negative = colgenie[colgenie.index.isin(all_negative)]
            dreamit_total_positive.append(coldreamit_positive.tolist())
            dreamit_total_negative.append(coldreamit_negative.tolist())
            dreamit_total.append([coldreamit_positive.tolist() + coldreamit_negative.tolist()])
            de_total_positive.append(colde_positive.tolist())
            de_total_negative.append(colde_negative.tolist())
            de_total.append([colde_positive.tolist() + colde_negative.tolist()])
            genie_total_positive.append(colgenie_positive.tolist())
            genie_total_negative.append(colgenie_negative.tolist())
            genie_total.append([colgenie_positive.tolist() + colgenie_negative.tolist()])
            data_details.append([col, tissue, species, len(markers2), len(all_negative)])

# plotting
# vals = [0, 1e-50, 1e-40, 1e-30, 1e-20, 1e-15, 1e-10, 1e-5, 1e-4, 1e-3, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1]
vals = [0, 1e-50, 1e-45, 1e-40, 1e-35, 1e-30, 1e-25, 1e-20, 1e-15, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3,
        0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05] + [i for i in np.arange(0.06, 1.01, 0.01)]
print('\n\n', 'vals', vals, '\n\n')
dreamit_TPR = []
dreamit_FPR = []
dreamit_precision = []
de_TPR = []
de_FPR = []
de_precision = []
genie_TPR = []
genie_FPR = []
genie_precision = []
for i in range(len(dreamit_total_positive)):
    dreamit_TPR_i = [0]
    dreamit_FPR_i = [0]
    dreamit_precision_i = [1]
    de_TPR_i = [0]
    de_FPR_i = [0]
    de_precision_i = [1]
    genie_TPR_i = [0]
    genie_FPR_i = [0]
    genie_precision_i = [1]
    for j in vals:
        dreamit_TP = len([x for x in dreamit_total_positive[i] if x <= j])
        dreamit_P = len(dreamit_total_positive[i])
        dreamit_FP = len([x for x in dreamit_total_negative[i] if x <= j])
        dreamit_N = len(dreamit_total_negative[i])
        tpr = dreamit_TP / dreamit_P  # also equal to recall
        fpr = dreamit_FP / dreamit_N
        if dreamit_TP > 0:
            precision = dreamit_TP / (dreamit_TP + dreamit_FP)
        else:
            precision = 1
        dreamit_TPR_i.append(tpr)  # also equal to recall
        dreamit_FPR_i.append(fpr)
        dreamit_precision_i.append(precision)

        de_TP = len([x for x in de_total_positive[i] if x <= j])
        de_P = len(de_total_positive[i])
        de_FP = len([x for x in de_total_negative[i] if x <= j])
        de_N = len(de_total_negative[i])
        detpr = de_TP / de_P
        defpr = de_FP / de_N
        if de_TP > 0:
            deprecision = de_TP / (de_TP + de_FP)
        else:
            deprecision = 1
        de_TPR_i.append(detpr)
        de_FPR_i.append(defpr)
        de_precision_i.append(deprecision)

        genie_TP = len([x for x in genie_total_positive[i] if x <= j])
        genie_P = len(genie_total_positive[i])
        genie_FP = len([x for x in genie_total_negative[i] if x <= j])
        genie_N = len(genie_total_negative[i])
        genietpr = genie_TP / genie_P
        geniefpr = genie_FP / genie_N
        if genie_TP > 0:
            genieprecision = genie_TP / (genie_TP + genie_FP)
        else:
            genieprecision = 1
        genie_TPR_i.append(genietpr)
        genie_FPR_i.append(geniefpr)
        genie_precision_i.append(genieprecision)

    # sorting
    zipped_lists = zip(dreamit_TPR_i, dreamit_FPR_i, dreamit_precision_i)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    dreamit_TPR_i, dreamit_FPR_i, dreamit_precision_i = [list(tuple) for tuple in tuples]

    zipped_lists = zip(de_TPR_i, de_FPR_i, de_precision_i)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    de_TPR_i, de_FPR_i, de_precision_i = [list(tuple) for tuple in tuples]

    zipped_lists = zip(genie_TPR_i, genie_FPR_i, genie_precision_i)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    genie_TPR_i, genie_FPR_i, genie_precision_i = [list(tuple) for tuple in tuples]

    # plotting ROC
    # Compute AUC using (y,x)
    auc_dreamit_i = np.trapz(dreamit_TPR_i, dreamit_FPR_i)
    auc_de_i = np.trapz(de_TPR_i, de_FPR_i)
    auc_genie_i = np.trapz(genie_TPR_i, genie_FPR_i)
    data_details[i].append(auc_dreamit_i)
    data_details[i].append(auc_de_i)
    data_details[i].append(auc_genie_i)
    plt.plot(dreamit_FPR_i, dreamit_TPR_i, linestyle='solid', color='blue',
             label=label1 + ' (AUC=%.2f)' % (auc_dreamit_i))
    plt.plot(de_FPR_i, de_TPR_i, linestyle='solid', color='red', label=label2 + ' (AUC=%.2f)' % (auc_de_i))
    plt.plot(genie_FPR_i, genie_TPR_i, linestyle='solid', color='green', label=label3 + ' (AUC=%.2f)' % (auc_genie_i))
    plt.scatter(dreamit_FPR_i, dreamit_TPR_i, color='blue')
    plt.scatter(de_FPR_i, de_TPR_i, color='red')
    plt.scatter(genie_FPR_i, genie_TPR_i, color='green')
    plt.plot([0., 1.], [0., 1.], linestyle='--', color='k')
    plt.legend(loc='lower right', fontsize=8)
    plt.title('ID: ' + data_details[i][0] + ' Tissue: ' + data_details[i][1] + ' Species: ' + data_details[i][
        2] + ' Markers Tested: ' + str(data_details[i][3]) + ' Others Tested: ' + str(data_details[i][4]))
    plt.ylabel('TPR')
    plt.xlabel('FPR')
    plt.savefig(output + 'new_roc_' + str(i) + '.png', dpi=300, bbox_inches='tight')
    plt.show()

    # plotting Precision-Recall
    # Compute AUC using (y,x)
    auc_dreamit_i = np.trapz(dreamit_precision_i, dreamit_TPR_i)
    auc_de_i = np.trapz(de_precision_i, de_TPR_i)
    auc_genie_i = np.trapz(genie_precision_i, genie_TPR_i)
    data_details[i].append(auc_dreamit_i)
    data_details[i].append(auc_de_i)
    data_details[i].append(auc_genie_i)
    plt.plot(dreamit_TPR_i, dreamit_precision_i, linestyle='solid', color='blue',
             label=label1 + ' (AUC=%.2f)' % (auc_dreamit_i))
    plt.plot(de_TPR_i, de_precision_i, linestyle='solid', color='red', label=label2 + ' (AUC=%.2f)' % (auc_de_i))
    plt.plot(genie_TPR_i, genie_precision_i, linestyle='solid', color='green',
             label=label3 + ' (AUC=%.2f)' % (auc_genie_i))
    plt.scatter(dreamit_TPR_i, dreamit_precision_i, color='blue')
    plt.scatter(de_TPR_i, de_precision_i, color='red')
    plt.scatter(genie_TPR_i, genie_precision_i, color='green')
    # plt.plot([0., 1.], [1., 0.], linestyle='--', color='k')
    plt.legend(loc='lower right', fontsize=8)
    plt.title('ID: ' + data_details[i][0] + ' Tissue: ' + data_details[i][1] + ' Species: ' + data_details[i][
        2] + ' Markers Tested: ' + str(data_details[i][3]) + ' Others Tested: ' + str(data_details[i][4]))
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    plt.savefig(output + 'new_pr_' + str(i) + '.png', dpi=300, bbox_inches='tight')
    plt.show()

    dreamit_TPR.append(dreamit_TPR_i)
    dreamit_FPR.append(dreamit_FPR_i)
    dreamit_precision.append(dreamit_precision_i)
    de_TPR.append(de_TPR_i)
    de_FPR.append(de_FPR_i)
    de_precision.append(de_precision_i)
    genie_TPR.append(de_TPR_i)
    genie_FPR.append(de_FPR_i)
    genie_precision.append(de_precision_i)

details_df = pd.DataFrame(data_details,
                          columns=['Branch', 'Tissue', 'Species', 'Factor Markers', 'Target Markers', 'ROC AUC DREAMIT2',
                                   'ROC AUC DE', 'ROC AUC GENIE3', 'PR AUC DREAMIT2', 'PR AUC DE', 'PR AUC GENIE3'])
print(details_df)
details_df.to_csv(output + 'summary.csv', index=False)

dreamit_TPR_mean = []
dreamit_FPR_mean = []
dreamit_precision_mean = []
de_TPR_mean = []
de_FPR_mean = []
de_precision_mean = []
genie_TPR_mean = []
genie_FPR_mean = []
genie_precision_mean = []
for i in range(len(dreamit_TPR[0])):
    dreamit_TPR_mean_i = np.mean([x[i] for x in dreamit_TPR])
    dreamit_FPR_mean_i = np.mean([x[i] for x in dreamit_FPR])
    dreamit_precision_mean_i = np.mean([x[i] for x in dreamit_precision])
    de_TPR_mean_i = np.mean([x[i] for x in de_TPR])
    de_FPR_mean_i = np.mean([x[i] for x in de_FPR])
    de_precision_mean_i = np.mean([x[i] for x in de_precision])
    genie_TPR_mean_i = np.mean([x[i] for x in genie_TPR])
    genie_FPR_mean_i = np.mean([x[i] for x in genie_FPR])
    genie_precision_mean_i = np.mean([x[i] for x in genie_precision])
    dreamit_TPR_mean.append(dreamit_TPR_mean_i)
    dreamit_FPR_mean.append(dreamit_FPR_mean_i)
    dreamit_precision_mean.append(dreamit_precision_mean_i)
    de_TPR_mean.append(de_TPR_mean_i)
    de_FPR_mean.append(de_TPR_mean_i)
    de_precision_mean.append(de_precision_mean_i)
    genie_TPR_mean.append(genie_TPR_mean_i)
    genie_FPR_mean.append(genie_TPR_mean_i)
    genie_precision_mean.append(genie_precision_mean_i)

# ROC plot
# Compute AUC using (y,x)
auc_dreamit = np.trapz(dreamit_TPR_mean, dreamit_FPR_mean)
auc_de = np.trapz(de_TPR_mean, de_FPR_mean)
auc_genie = np.trapz(genie_TPR_mean, genie_FPR_mean)
plt.plot(dreamit_FPR_mean, dreamit_TPR_mean, linestyle='solid', color='blue',
         label=label1 + ' (AUC=%.2f)' % (auc_dreamit))
plt.plot(de_FPR_mean, de_TPR_mean, linestyle='solid', color='red', label=label2 + ' (AUC=%.2f)' % (auc_de))
plt.plot(genie_FPR_mean, genie_TPR_mean, linestyle='solid', color='green', label=label3 + ' (AUC=%.2f)' % (auc_genie))
plt.scatter(dreamit_FPR_mean, dreamit_TPR_mean, color='blue')
plt.scatter(de_FPR_mean, de_TPR_mean, color='red')
plt.scatter(genie_FPR_mean, genie_TPR_mean, color='green')
plt.plot([0., 1.], [0., 1.], linestyle='--', color='k')
plt.legend(loc='lower right', fontsize=8)
plt.title('Mean TPR and FPR Observations for ' + label1 + ', ' + label2 + ', and ' + label3)
plt.ylabel('TPR')
plt.xlabel('FPR')
plt.savefig(output + 'new_roc_allmetrics_mean_nomouse.png', dpi=300)
plt.show()

# Precision-Recall plot
# Compute AUC using (y,x)
auc_dreamit = np.trapz(dreamit_precision_mean, dreamit_TPR_mean)
auc_de = np.trapz(de_precision_mean, de_TPR_mean)
auc_genie = np.trapz(genie_precision_mean, genie_TPR_mean)
plt.plot(dreamit_TPR_mean, dreamit_precision_mean, linestyle='solid', color='blue',
         label=label1 + ' (AUC=%.2f)' % (auc_dreamit))
plt.plot(de_TPR_mean, de_precision_mean, linestyle='solid', color='red', label=label2 + ' (AUC=%.2f)' % (auc_de))
plt.plot(genie_TPR_mean, genie_precision_mean, linestyle='solid', color='green',
         label=label3 + ' (AUC=%.2f)' % (auc_genie))
plt.scatter(dreamit_TPR_mean, dreamit_precision_mean, color='blue')
plt.scatter(de_TPR_mean, de_precision_mean, color='red')
plt.scatter(genie_TPR_mean, genie_precision_mean, color='green')
# plt.plot([0., 1.], [1., 0.], linestyle='--', color='k')
plt.legend(loc='lower right', fontsize=8)
plt.title('Mean Precision-Recall Observations for ' + label1 + ', ' + label2 + ', and ' + label3)
plt.ylabel('Precision')
plt.xlabel('Recall')
plt.savefig(output + 'new_pr_allmetrics_mean_nomouse.png', dpi=300)
plt.show()

# aggregate all
# flat_list = [item for sublist in t for item in sublist]
dreamit_agg_positive = [i for sub in dreamit_total_positive for i in sub]
dreamit_agg_negative = [i for sub in dreamit_total_negative for i in sub]
de_agg_positive = [i for sub in de_total_positive for i in sub]
de_agg_negative = [i for sub in de_total_negative for i in sub]
genie_agg_positive = [i for sub in genie_total_positive for i in sub]
genie_agg_negative = [i for sub in genie_total_negative for i in sub]

# plotting
# vals = [0, 1e-5, 1e-4, 1e-3, 0.01, 0.02, 0.03, 0.04, 0.05] + [i for i in np.arange(0.06, 1.01, 0.01)]
dreamit_TPR = []
dreamit_FPR = []
dreamit_precision = []
de_TPR = []
de_FPR = []
de_precision = []
genie_TPR = []
genie_FPR = []
genie_precision = []
for j in vals:
    dreamit_TP = len([x for x in dreamit_agg_positive if x <= j])
    dreamit_P = len(dreamit_agg_positive)
    dreamit_FP = len([x for x in dreamit_agg_negative if x <= j])
    dreamit_N = len(dreamit_agg_negative)
    tpr = dreamit_TP / dreamit_P
    fpr = dreamit_FP / dreamit_N
    if dreamit_TP > 0:
        precision = dreamit_TP / (dreamit_TP + dreamit_FP)
    else:
        precision = 1
    dreamit_TPR.append(tpr)
    dreamit_FPR.append(fpr)
    dreamit_precision.append(precision)

    de_TP = len([x for x in de_agg_positive if x <= j])
    de_P = len(de_agg_positive)
    de_FP = len([x for x in de_agg_negative if x <= j])
    de_N = len(de_agg_negative)
    detpr = de_TP / de_P
    defpr = de_FP / de_N
    if de_TP > 0:
        deprecision = de_TP / (de_TP + de_FP)
    else:
        deprecision = 1
    de_TPR.append(detpr)
    de_FPR.append(defpr)
    de_precision.append(deprecision)

    genie_TP = len([x for x in genie_agg_positive if x <= j])
    genie_P = len(genie_agg_positive)
    genie_FP = len([x for x in genie_agg_negative if x <= j])
    genie_N = len(genie_agg_negative)
    genietpr = genie_TP / genie_P
    geniefpr = genie_FP / genie_N
    if genie_TP > 0:
        genieprecision = genie_TP / (genie_TP + genie_FP)
    else:
        genieprecision = 1
    genie_TPR.append(genietpr)
    genie_FPR.append(geniefpr)
    genie_precision.append(genieprecision)

# ROC plot
# Compute AUC using (y,x)
roc_auc_dreamit = np.trapz(dreamit_TPR, dreamit_FPR)
roc_auc_de = np.trapz(de_TPR, de_FPR)
roc_auc_genie = np.trapz(genie_TPR, genie_FPR)
plt.plot(dreamit_FPR, dreamit_TPR, linestyle='solid', color='blue', label=label1 + ' (AUC=%.2f)' % (roc_auc_dreamit))
plt.plot(de_FPR, de_TPR, linestyle='solid', color='red', label=label2 + ' (AUC=%.2f)' % (roc_auc_de))
plt.plot(genie_FPR, genie_TPR, linestyle='solid', color='green', label=label3 + ' (AUC=%.2f)' % (roc_auc_genie))
plt.scatter(dreamit_FPR, dreamit_TPR, color='blue')
plt.scatter(de_FPR, de_TPR, color='red')
plt.scatter(genie_FPR, genie_TPR, color='green')
plt.plot([0., 1.], [0., 1.], linestyle='--', color='k')
plt.legend(loc='lower right', fontsize=8)
plt.title('Aggregate TPR and FPR Observations for ' + label1 + ', ' + label2 + ', and ' + label3)
plt.ylabel('TPR')
plt.xlabel('FPR')
plt.savefig(output + 'new_roc_allmetrics_agg_nomouse.png', dpi=300)
plt.show()

# Precision-Recall plot
# Compute AUC using (y,x)
pr_auc_dreamit = np.trapz(dreamit_precision, dreamit_TPR)
pr_auc_de = np.trapz(de_precision, de_TPR)
pr_auc_genie = np.trapz(genie_precision, genie_TPR)
plt.plot(dreamit_TPR, dreamit_precision, linestyle='solid', color='blue',
         label=label1 + ' (AUC=%.2f)' % (pr_auc_dreamit))
plt.plot(de_TPR, de_precision, linestyle='solid', color='red', label=label2 + ' (AUC=%.2f)' % (pr_auc_de))
plt.plot(genie_TPR, genie_precision, linestyle='solid', color='green', label=label3 + ' (AUC=%.2f)' % (pr_auc_genie))
plt.scatter(dreamit_TPR, dreamit_precision, color='blue')
plt.scatter(de_TPR, de_precision, color='red')
plt.scatter(genie_TPR, genie_precision, color='green')
# plt.plot([0., 1.], [1., 0.], linestyle='--', color='k')
plt.legend(loc='lower right', fontsize=8)
plt.title('Aggregate Precision-Recall Observations for ' + label1 + ', ' + label2 + ', and ' + label3)
plt.ylabel('Precision')
plt.xlabel('Recall')
plt.savefig(output + 'new_pr_allmetrics_agg_nomouse.png', dpi=300)
plt.show()

print(label1 + 'TPR: ', dreamit_TPR)
print(label1 + 'FPR: ', dreamit_FPR)
print(label1 + 'precision:', dreamit_precision)

print(label2 + 'TPR: ', de_TPR)
print(label2 + 'FPR: ', de_FPR)
print(label2 + 'precision:', de_precision)

print(label3 + 'TPR: ', genie_TPR)
print(label3 + 'FPR: ', genie_FPR)
print(label3 + 'precision:', genie_precision)
