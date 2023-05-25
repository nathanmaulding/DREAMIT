import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import os

# for roc_data2.py
dreamit_df = pd.read_csv('/Users/nathandmaulding/Desktop/specificity_results/v20/fdr/dreamit_fdr_results_v20.csv', index_col=0) # old data: dreamit_fdr_results3.csv
de_df = pd.read_csv('/Users/nathandmaulding/Desktop/specificity_results/v20/fdr/de_fdr_results_v20_raw.csv', index_col=0) # old data: de_fdr_results3.csv
genie_df = pd.read_csv('/Users/nathandmaulding/Desktop/specificity_results/v20/fdr/genie_fdr_results.csv', index_col=0)
output = '/Users/nathandmaulding/Desktop/specificity_results/v20/specificity_universe/'
label1 = 'DREAMIT2'
label2 = 'DE raw'
label3 = 'GENIE3'

# unique factors
Testis = ['TFAP2C', 'PRDM1', 'PRDM14', 'EPHA2', 'PLXDC2', 'RND3', 'TIMP3', 'NR2F1', 'NR2F2', 'DMRT1', 'ADGRA3', 'ID4', 'KIT', 'UTF1', 'TAF4B']
Brain = ['GAP43', 'IGF1', 'H19', 'CTCF', 'MLANA', 'PAX3', 'GLI3', 'MAP2', 'HPRT1', 'MET', 'SLC7A2', 'CDC25B', 'PLK1', 'BDNF', 'NGF', 'CDH5', 'FOXG1', 'GSK3B', 'KLF9', 'ILF2', 'HSPA8', 'ZNF423', 'REST', 'CDKN1C', 'THRAP3', 'PIK3CA', 'OPA1', 'FLI1', 'EWSR1', 'CD99', 'IDH1', 'DVL1', 'DVL2', 'DVL3', 'MMP9', 'PDGFRA', 'NFE2L2', 'NODAL', 'ATM', 'HINT1', 'CAST', 'FABP7', 'ADGRB1', 'ELOVL2', 'ARNT2', 'SALL2', 'PKNOX1', 'EIF4E2', 'DDX28', 'FAT1', 'BCL2L11', 'MIR486-1', 'NR4A2', 'ASCL1', 'LMX1B', 'TH', 'DDC', 'S1PR5', 'IL1B', 'SLC16A1', 'SOD1', 'KRIT1', 'PDCD10', 'MC4R', 'SOX11', 'SOX12', 'FGF8', 'CFH', 'SIRPB1', 'IGFBP5', 'VRK2', 'BAX', 'PAX4', 'ZKSCAN7', 'MAPT', 'NPAS4', 'LHX4', 'PTK2', 'SOX10', 'MYRF', 'PBX1', 'ZBTB24', 'MCM5', 'PCNA', 'CDK2', 'CDC25A']
Bone_marrow = ['CEBPA', 'CD7', 'FUT4', 'FLT3', 'BAALC', 'MAFB', 'EGFR', 'PTTG1', 'MAPK7', 'IL3RA', 'IKZF1', 'IKZF3', 'JAK1', 'JAK2', 'RMND1', 'RFX3', 'BCR', 'RUNX1T1', 'NEDD8', 'TAZ', 'LATS1', 'DLX2', 'HHEX', 'HLX', 'HMX1', 'NKX31', 'DLX5', 'DLX6', 'TET2', 'POU2AF1', 'EDIL3', 'NSMAF', 'TPSAB1', 'CD34', 'HOXB4', 'IL3', 'IL11', 'SREBF1', 'SOX5', 'SOX6', 'E2F2', 'LDB2', 'MECOM', 'EBF1', 'PAX5', 'BMP2', 'ALPL', 'BGLAP', 'IL17A', 'MIR15B', 'SMURF1', 'ACAN', 'COL2A1', 'EGR3', 'DSPP', 'IBSP', 'ATXN1', 'ITGB1', 'ACE', 'SNW1', 'SHH', 'FUBP1', 'TAL1', 'KDM1A', 'GFI1B', 'IRF8', 'MOV10', 'AGO3', 'HDAC1', 'RBBP5', 'WASF2', 'SDC1', 'IL2', 'TRPM7', 'MAGT1', 'WNT11', 'KLF7', 'OSTF1', 'IRF7', 'MIR223', 'IFNB1', 'IRF5', 'IRF9', 'WNT1', 'DNMT3A', 'BMP6', 'SMAD1', 'SMAD5', 'SMAD9', 'ELANE', 'MPO', 'GFI1', 'CEBPD', 'CEBPE']
Liver = ['FOXA1', 'HNF4A', 'HNF1B', 'CDKN2A', 'ID1', 'TTF1', 'CYP3A4', 'UGT1A1', 'ABCC2', 'FOXA2', 'CYP7A1', 'AXIN1', 'AKR1C2', 'AKR1C3', 'EPCAM', 'SLC5A8', 'NR1I2', 'ALBUMINFAMILY', 'CXCL8', 'GPX4', 'RUNX3', 'VWCE', 'MIR148A', 'AFP', 'CD4', 'KRT19', 'FN1', 'TWIST2', 'ABCC1', 'GLI2', 'BVES', 'SLCO1B3', 'ZEB2', 'UGT1A3', 'UGT2B17', 'LHX2', 'ABCG1', 'ABCG2', 'YAP1', 'TEAD4', 'NR1H4', 'DGKB', 'BPTF', 'CDX2', 'MANBA', 'UBE2M', 'GPER1', 'P2RY11', 'ARG1', 'PKLR', 'PCSK9', 'PNPLA3', 'PPIG', 'RBM8A', 'HMGB1', 'CYP2D6', 'CYP8B1', 'FOXN3', 'PKM', 'PTEN', 'CD274', 'MIR214', 'SLC7A5', 'ZBTB7C', 'PCK1', 'G6PC', 'SLC46A3', 'USF1', 'SLC25A13', 'TBP', 'SP3', 'ATP2A3', 'CDNF', 'CRELD2', 'DNAJB11', 'DTL', 'GINS2', 'MANF', 'PDIA4', 'PDIA6', 'VCP', 'PHLPP1', 'SMARCA2', 'GATA3', 'MIR29A', 'CD36', 'DDX17', 'DAB2IP', 'HAVCR2', 'CTLA4', 'BATF', 'EOMES', 'XRCC5', 'PCK2', 'RREB1', 'ETS1', 'IRF1', 'IRF2', 'E2F7', 'E2F8', 'MYCN', 'IRS1', 'IRS2', 'KMT5A', 'TFCP2', 'ZHX2', 'KDM2A', 'GNAS', 'CD28', 'SYP', 'KRT7', 'GFAP', 'ONECUT1', 'ZP4', 'MAZ', 'CTNNBIP1', 'LRAT', 'SMAD3', 'MLXIPL', 'ADCY3', 'GTF3A', 'TNFRSF10B', 'HSPA5', 'NCAM1', 'CXCR4', 'LGR5', 'NR3C1', 'NQO1', 'HGF', 'ERAL1', 'IL12RB2', 'TIMP1', 'FOXA3', 'EPAS1', 'NR1H3', 'NR1H2', 'TET3', 'MAVS', 'TFE3', 'MAF', 'NFE2L3', 'FOS', 'ARNT', 'ZC3H12A', 'TFEC', 'MEIS2']
Retina = ['E2F4', 'RBL2', 'TJP1', 'GDF10', 'RAX', 'VSX2', 'PAX6', 'CPNE4', 'CPNE5', 'CPNE6', 'PITX2', 'ATOH7', 'OTX2']
Heart = ['HAND1', 'HAND2', 'TBX5', 'ESRRG', 'MESP1', 'NKX2-5', 'ZFPM2', 'APEX1', 'DECR1', 'ALOX15', 'TTN', 'TCF21', 'WT1', 'TBX18', 'ACTC1', 'MCL1', 'TBX20', 'NOTCH1', 'NOTCH2', 'PECAM1', 'KDR', 'FGFR1', 'VWF', 'HEY2', 'KDM5A', 'MAPK8', 'ESR1']
Universe = Testis + Brain + Bone_marrow + Liver + Retina + Heart
UniDict = {'Testis': Testis, 'Brain': Brain, 'Bone marrow': Bone_marrow, 'Liver': Liver, 'Retina': Retina, 'Heart': Heart}

dreamit_df.columns = [i.replace('_v20', '') for i in dreamit_df.columns.tolist()]
de_df.columns = [i.replace('_v20', '') for i in de_df.columns.tolist()]

tissue_dict = {'EBI3': 'Testis', 'EBI4': 'Brain', 'EBI5': 'Bone marrow', 'EBI6': 'Liver', 'EBI7': 'Embryo', 'EBI8': 'Embryo', 'EBI9': 'Retina', 'EBI10': 'Brain', 'friedman': 'Heart', 'paul': 'Embryo'}

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
        if each in col:  #and 'EBI7' not in col and 'EBI3' not in col and 'EBI9' not in col:  # added and 'paul' not in col
            tissue = tissue_dict[each]
            factor_markers = UniDict[tissue]
            print('Tissue: ', tissue, factor_markers)
            if 'EBI' in col or 'friedman' in col:
                species = 'Human'
            else:
                species = 'Mouse'
            # for DE
            de_df = de_df[de_df.index.isin(Universe)]
            sub_de_df = de_df[de_df.index.isin(factor_markers)]
            if len(sub_de_df[col].dropna().index) > 0:
                markers = sub_de_df[col].dropna().index.tolist()
                not_sub = de_df[col].drop(markers).dropna()
                de_total.append(de_df[col].dropna().tolist())
                de_total_positive.append(sub_de_df[col].dropna().tolist())
                de_total_negative.append(not_sub.tolist())
                de_total_factors.append(sub_de_df.index.tolist())
            # for DREAMIT2
            dreamit_df = dreamit_df[dreamit_df.index.isin(Universe)]
            sub_dreamit_df = dreamit_df[dreamit_df.index.isin(factor_markers)]
            if len(sub_dreamit_df[col].dropna().index) > 0:
                markers = sub_dreamit_df[col].dropna().index.tolist()
                not_sub = dreamit_df[col].drop(markers).dropna()
                dreamit_total.append(dreamit_df[col].dropna().tolist())
                dreamit_total_positive.append(sub_dreamit_df[col].dropna().tolist())
                dreamit_total_negative.append(not_sub.tolist())
                dreamit_total_factors.append(sub_dreamit_df.index.tolist())
                data_details.append([col, tissue, species, len(sub_dreamit_df[col].dropna().tolist()), len(not_sub.tolist())])
            # for GENIE3
            genie_df = genie_df[genie_df.index.isin(Universe)]
            sub_genie_df = genie_df[genie_df.index.isin(factor_markers)]
            # if equal factor numbers between methods necessary use line below
            #sub_genie_df = sub_genie_df[sub_genie_df.index.isin(sub_dreamit_df.index.tolist())]
            if len(sub_genie_df[col].dropna().index) > 0:
                markers = sub_genie_df[col].dropna().index.tolist()
                not_sub = genie_df[col].drop(markers).dropna()
                genie_total.append(genie_df[col].dropna().tolist())
                genie_total_positive.append(sub_genie_df[col].dropna().tolist())
                genie_total_negative.append(not_sub.tolist())
                genie_total_factors.append(sub_genie_df.index.tolist())

print(len(data_details), data_details)
print(len(dreamit_total), dreamit_total)
print(len(dreamit_total_positive), dreamit_total_positive)
print(len(dreamit_total_negative), dreamit_total_negative)
print(len(de_total), de_total)
print(len(de_total_positive), de_total_positive)
print(len(de_total_negative), de_total_negative)
print(len(genie_total), genie_total)
print(len(genie_total_positive), genie_total_positive)
print(len(genie_total_negative), genie_total_negative)

defax = []
dreamitfax = []
geniefax = []
for each in de_total_factors:
    for i in each:
        if i not in defax:
            defax.append(i)
for each in dreamit_total_factors:
    for i in each:
        if i not in dreamitfax:
            dreamitfax.append(i)
for each in genie_total_factors:
    for i in each:
        if i not in geniefax: #and i in dreamitfax: # has to be in dreamit to keep fax #s equal/comparable?
            geniefax.append(i)

print(len(defax), defax)
print(len(dreamitfax), dreamitfax)
print(len(geniefax), geniefax)


# plotting
#vals = [0, 1e-50, 1e-40, 1e-30, 1e-20, 1e-15, 1e-10, 1e-5, 1e-4, 1e-3, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1]
vals = [0, 1e-50, 1e-45, 1e-40, 1e-35, 1e-30, 1e-25, 1e-20, 1e-15, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05] + [i for i in np.arange(0.06, 1.01, 0.01)]
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
        tpr = dreamit_TP / dreamit_P # also equal to recall
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
    plt.plot(dreamit_FPR_i, dreamit_TPR_i, linestyle='solid', color='blue', label=label1+ ' (AUC=%.2f)' % (auc_dreamit_i))
    plt.plot(de_FPR_i, de_TPR_i, linestyle='solid', color='red', label=label2+' (AUC=%.2f)' % (auc_de_i))
    plt.plot(genie_FPR_i, genie_TPR_i, linestyle='solid', color='green', label=label3 + ' (AUC=%.2f)' % (auc_genie_i))
    plt.scatter(dreamit_FPR_i, dreamit_TPR_i, color='blue')
    plt.scatter(de_FPR_i, de_TPR_i, color='red')
    plt.scatter(genie_FPR_i, genie_TPR_i, color='green')
    plt.plot([0., 1.], [0., 1.], linestyle='--', color='k')
    plt.legend(loc='lower right', fontsize=8)
    plt.title('ID: ' + data_details[i][0] + ' Tissue: ' + data_details[i][1] + ' Species: ' + data_details[i][2] + ' Markers Tested: ' + str(data_details[i][3]) + ' Others Tested: ' + str(data_details[i][4]))
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
    plt.plot(dreamit_TPR_i, dreamit_precision_i, linestyle='solid', color='blue', label=label1+' (AUC=%.2f)' % (auc_dreamit_i))
    plt.plot(de_TPR_i, de_precision_i, linestyle='solid', color='red', label=label2+' (AUC=%.2f)' % (auc_de_i))
    plt.plot(genie_TPR_i, genie_precision_i, linestyle='solid', color='green', label=label3 + ' (AUC=%.2f)' % (auc_genie_i))
    plt.scatter(dreamit_TPR_i, dreamit_precision_i, color='blue')
    plt.scatter(de_TPR_i, de_precision_i, color='red')
    plt.scatter(genie_TPR_i, genie_precision_i, color='green')
    #plt.plot([0., 1.], [1., 0.], linestyle='--', color='k')
    plt.legend(loc='lower right', fontsize=8)
    plt.title('ID: ' + data_details[i][0] + ' Tissue: ' + data_details[i][1] + ' Species: ' + data_details[i][2] + ' Markers Tested: ' + str(data_details[i][3]) + ' Others Tested: ' + str(data_details[i][4]))
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

details_df = pd.DataFrame(data_details, columns=['Branch', 'Tissue', 'Species', 'Factor Markers', 'Target Markers', 'ROC AUC DREAMIT2', 'ROC AUC DE', 'ROC AUC GENIE3', 'PR AUC DREAMIT2', 'PR AUC DE', 'PR AUC GENIE3'])
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
plt.plot(dreamit_FPR_mean, dreamit_TPR_mean, linestyle='solid', color='blue', label=label1+' (AUC=%.2f)' % (auc_dreamit))
plt.plot(de_FPR_mean, de_TPR_mean, linestyle='solid', color='red', label=label2+' (AUC=%.2f)' % (auc_de))
plt.plot(genie_FPR_mean, genie_TPR_mean, linestyle='solid', color='green', label=label3+' (AUC=%.2f)' % (auc_genie))
plt.scatter(dreamit_FPR_mean, dreamit_TPR_mean, color='blue')
plt.scatter(de_FPR_mean, de_TPR_mean, color='red')
plt.scatter(genie_FPR_mean, genie_TPR_mean, color='green')
plt.plot([0., 1.], [0., 1.], linestyle='--', color='k')
plt.legend(loc='lower right', fontsize=8)
plt.title('Mean TPR and FPR Observations for '+label1+', '+label2+', and '+label3)
plt.ylabel('TPR')
plt.xlabel('FPR')
plt.savefig(output + 'new_roc_allmetrics_mean_nomouse.png', dpi=300)
plt.show()

# Precision-Recall plot
# Compute AUC using (y,x)
auc_dreamit = np.trapz(dreamit_precision_mean, dreamit_TPR_mean)
auc_de = np.trapz(de_precision_mean, de_TPR_mean)
auc_genie = np.trapz(genie_precision_mean, genie_TPR_mean)
plt.plot(dreamit_TPR_mean, dreamit_precision_mean, linestyle='solid', color='blue', label=label1+' (AUC=%.2f)' % (auc_dreamit))
plt.plot(de_TPR_mean, de_precision_mean, linestyle='solid', color='red', label=label2+' (AUC=%.2f)' % (auc_de))
plt.plot(genie_TPR_mean, genie_precision_mean, linestyle='solid', color='green', label=label3+' (AUC=%.2f)' % (auc_genie))
plt.scatter(dreamit_TPR_mean, dreamit_precision_mean, color='blue')
plt.scatter(de_TPR_mean, de_precision_mean, color='red')
plt.scatter(genie_TPR_mean, genie_precision_mean, color='green')
#plt.plot([0., 1.], [1., 0.], linestyle='--', color='k')
plt.legend(loc='lower right', fontsize=8)
plt.title('Mean Precision-Recall Observations for '+label1+', '+label2+', and '+label3)
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
#vals = [0, 1e-5, 1e-4, 1e-3, 0.01, 0.02, 0.03, 0.04, 0.05] + [i for i in np.arange(0.06, 1.01, 0.01)]
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
plt.plot(dreamit_FPR, dreamit_TPR, linestyle='solid', color='blue', label=label1+' (AUC=%.2f)' % (roc_auc_dreamit))
plt.plot(de_FPR, de_TPR, linestyle='solid', color='red', label=label2+' (AUC=%.2f)' % (roc_auc_de))
plt.plot(genie_FPR, genie_TPR, linestyle='solid', color='green', label=label3+' (AUC=%.2f)' % (roc_auc_genie))
plt.scatter(dreamit_FPR, dreamit_TPR, color='blue')
plt.scatter(de_FPR, de_TPR, color='red')
plt.scatter(genie_FPR, genie_TPR, color='green')
plt.plot([0., 1.], [0., 1.], linestyle='--', color='k')
plt.legend(loc='lower right', fontsize=8)
plt.title('Aggregate TPR and FPR Observations for '+label1+', '+label2+', and '+label3)
plt.ylabel('TPR')
plt.xlabel('FPR')
plt.savefig(output + 'new_roc_allmetrics_agg_nomouse.png', dpi=300)
plt.show()


# Precision-Recall plot
# Compute AUC using (y,x)
pr_auc_dreamit = np.trapz(dreamit_precision, dreamit_TPR)
pr_auc_de = np.trapz(de_precision, de_TPR)
pr_auc_genie = np.trapz(genie_precision, genie_TPR)
plt.plot(dreamit_TPR, dreamit_precision, linestyle='solid', color='blue', label=label1+' (AUC=%.2f)' % (pr_auc_dreamit))
plt.plot(de_TPR, de_precision, linestyle='solid', color='red', label=label2+' (AUC=%.2f)' % (pr_auc_de))
plt.plot(genie_TPR, genie_precision, linestyle='solid', color='green', label=label3+' (AUC=%.2f)' % (pr_auc_genie))
plt.scatter(dreamit_TPR, dreamit_precision, color='blue')
plt.scatter(de_TPR, de_precision, color='red')
plt.scatter(genie_TPR, genie_precision, color='green')
#plt.plot([0., 1.], [1., 0.], linestyle='--', color='k')
plt.legend(loc='lower right', fontsize=8)
plt.title('Aggregate Precision-Recall Observations for '+label1+', '+label2+', and '+label3)
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