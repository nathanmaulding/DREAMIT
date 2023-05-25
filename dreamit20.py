# ptime smoothing
########### CommandLine ############################
import pandas as pd


class CommandLine():

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.

        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        # branch, number_of_bins, bin_by_var, even_distribution, tf_dict
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )
        self.parser.add_argument('-json_dictionary',
                                 help="input a JSON file with branching information (created by NM_slingshot_v2.R)")
        self.parser.add_argument('-tf_dict',
                                 help="input must be a factor to all targets file with scores (could be edited)")
        self.parser.add_argument('-threads', type=int, default=1,
                                 help='input the number of threads to use. Default = 1')
        self.parser.add_argument('-number_of_bins', type=int, default=None,
                                 help='Input option for the number of bins to use (default = 20)')
        self.parser.add_argument('-bin_by_var', type=str, default="psuedotime",
                                 help='Option to bin by any string. Not developed, leave as default="psuedotime"')
        self.parser.add_argument('-even_distribution', type=bool, default=False,
                                 help='Option for creating bins with an equal distribution of cells in each bin. Not developed, (default = False)')
        self.parser.add_argument('-outdir', type=str, default='./dreamit_analysis')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


##################################################

class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''

    def __init__(self, msg):
        self.msg = msg


#####################################################

import json
import numpy as np
import scipy.stats as stats
import math
import matplotlib.pyplot as plt
from scipy import interpolate
import random
from scipy.stats import ks_2samp
import operator
from sklearn import metrics
from dtw import dtw, accelerated_dtw
import statsmodels.stats.multitest
from multiprocessing import Pool
from itertools import repeat
import sys
import os
import time


class binData:
    '''This function flows from init -> smoother -> binningData[bin] -> getData -> avgExpression
    -> lookupCells -> make avgGeneDict binningData[bin] -> next bin -> smoother

    Written by: Maulding
    '''

    def __init__(self, branch, outdir, number_of_bins, bin_by_var, even_distribution, tf_dict, data_dict):
        self.branch = branch
        self.outdir = outdir
        self.number_of_bins = number_of_bins
        self.bin_by_var = bin_by_var
        self.even_distribution = even_distribution
        self.data_dict = data_dict
        self.tf_dict = tf_dict
        # self.smoother()

    def getData(self):
        self.cells = self.data_dict[self.branch]["cell_id"]  # get the index of each cell
        self.cell_var = [{} for _ in self.cells]  # should be cell: var dictionary
        self.var_cell = [{} for _ in self.var]  # should be a var: cell dictionary
        for index, time in enumerate(self.var):  # iter through the list of vars
            self.var[index] = float(time)  # change vars to float
            self.cell_var[index][self.cells[index]] = time  # dict of cells (key) and time (vals)
            self.var_cell[index][self.var[index]] = self.cells[index]  # dict of time (keys) and cells (vals)

    def lookupCells(self, cells_in_bin, gene):
        # requires the list of cells in the bin and the gene as arguments
        sumGene = 0
        for cell_info in cells_in_bin:
            # use the index of the cell in question to get the gene expression
            sumGene += float(self.data_dict[self.branch][gene][cell_info[0]])
        return sumGene / len(cells_in_bin)

    def lookupCells2(self, cells_in_bin, gene, rawdf):
        # requires the list of cells in the bin and the gene as arguments
        sumGene = 0
        # print(rawdf)
        # print(gene, cells_in_bin)
        # print(rawdf[gene])
        # print(rawdf[gene].tolist())
        sys.stdout.flush()
        for cell_info in cells_in_bin:
            # use the index of the cell in question to get the gene expression
            if cell_info[0] not in rawdf[gene]:
                print('not in')
                print(rawdf[gene])
                print(rawdf[gene].tolist())
                print(rawdf[gene].index.tolist())
                print(cell_info[0])
                sys.stdout.flush()
            sumGene += float(rawdf[gene][cell_info[0]])
        return sumGene / len(cells_in_bin)

    def avgExpression(self, bin_index, cells_in_this_bin):
        if len(cells_in_this_bin) > 0:
            self.avgGeneBinDict[bin_index] = {}
            for each_gene in list(self.data_dict[self.branch].keys())[9:]:  # should be only genes
                self.avgGeneBinDict[bin_index][each_gene] = self.lookupCells(cells_in_this_bin, each_gene)
            for i, rfact in enumerate(self.randfactors):
                self.avgGeneBinDict[bin_index]['Random' + str(i + 1)] = self.lookupCells(cells_in_this_bin, rfact)
        # should return average expression of the gene expression of the cells in that bin
        # shouldn't even need to return it #return self.avgGeneBinDict

    def binningData(self):
        # used for smoothing by equal ptime intervals
        self.getData()
        self.avgGeneBinDict = {}
        bin_range = self.var_max - self.var_min  # get min and max time on branch
        bin_increment = bin_range / self.number_of_bins  # create an increment to move by
        bin_min = self.var_min
        self.bins = [self.var_min]  # create bins based on minimum increments
        self.bin_dict = []
        random.seed(42)
        self.randfactors = random.choices(
            [k for k in self.tf_dict.keys() if k in list(self.data_dict[self.branch].keys())[9:]], k=3)
        for i in range(self.number_of_bins):
            self.bin_dict.append({})
            bin_min += bin_increment
            self.bins.append(bin_min)  # span bins by the the previous bins min + the increment
        for index, timepoint in enumerate(self.bins):  # for the index (which is the minimum time) of the bin
            cells_in_this_bin = []
            if index >= len(self.bins) - 1:  # don't go over range
                pass
            else:
                for i, each in enumerate(self.data_dict[self.branch][self.bin_by_var]):
                    # if the cell's ptime is >= the bin, but < bin + 1
                    # then it is within the range/increment and is part of that bin
                    if float(each) >= timepoint and float(each) < self.bins[index + 1]:
                        cells_in_this_bin.append((i, each))
                # pass the bin and cells in the bin for averaging
                self.avgExpression(index, cells_in_this_bin)
        self.avgGeneBinDict['psuedotime_med'] = []
        for i in range(len(self.bins) - 1):
            self.avgGeneBinDict['psuedotime_med'].append(np.mean([self.bins[i], self.bins[i + 1]]))
        return self.avgGeneBinDict

    def binningData2(self):
        # used for smoothing by equal cell num
        self.getData()
        self.avgGeneBinDict = {}
        random.seed(42)
        self.randfactors = random.choices(
            [k for k in self.tf_dict.keys() if k in list(self.data_dict[self.branch].keys())[9:]], k=3)
        cellnum = len(self.var)
        d = np.floor(cellnum / self.number_of_bins)
        rem = cellnum - (d * self.number_of_bins)
        bin_content = [int(d + 1) if i < rem else int(d) for i in range(self.number_of_bins)]
        self.avgGeneBinDict = {bin: {} for bin in range(self.number_of_bins)}
        for each_gene in list(self.data_dict[self.branch].keys())[9:]:
            rawexp = [float(i) for i in self.data_dict[self.branch][each_gene]]
            rawptime2, rawexp2 = zip(*sorted(zip(self.var, rawexp)))
            iter = 0
            last = 0
            while iter < len(bin_content):  # number of bins
                self.avgGeneBinDict[iter][each_gene] = np.mean(rawexp2[last:last + bin_content[iter]])
                last += bin_content[iter]
                iter += 1
        for j, rfact in enumerate(self.randfactors):
            rawexp = [float(i) for i in self.data_dict[self.branch][rfact]]
            cellids = self.data_dict[self.branch]['cell_id']
            rawptime2, rawexp2, cellids2 = zip(*sorted(zip(self.var, rawexp, cellids)))
            iter = 0
            last = 0
            while iter < len(bin_content):  # number of bins
                self.avgGeneBinDict[iter]['Random' + str(j + 1)] = np.mean(rawexp2[last:last + bin_content[iter]])
                last += bin_content[iter]
                iter += 1
            # write code to export which cells are in each bin here
            # and averaged ptime as self object for use in last lines of smoother function
            if j == 1:  # run just once
                iterj = 0
                lastj = 0
                # print('bin content', bin_content, sum(bin_content), len(self.data_dict[self.branch]['cell_id']))
                # print('first slice', rawptime2[lastj:lastj + bin_content[iterj]], iterj, lastj, bin_content[iterj])
                while iterj < len(bin_content):  # number of bins
                    self.avgGeneBinDict[iterj]['psuedotime_mean'] = np.mean(rawptime2[lastj:lastj + bin_content[iterj]])
                    self.avgGeneBinDict[iterj]['psuedotime_med'] = np.median(
                        rawptime2[lastj:lastj + bin_content[iterj]])
                    self.avgGeneBinDict[iterj]['cell_ids'] = cellids2[lastj:lastj + bin_content[iterj]]
                    lastj += bin_content[iterj]
                    iterj += 1
        return self.avgGeneBinDict

    def makeBinList(self, bin):
        bin_range = self.var_max - self.var_min  # get min and max time on branch
        bin_increment = bin_range / bin  # create an increment to move by
        bin_min = self.var_min
        bins = [bin_min]
        for i in range(bin):
            bin_min += bin_increment
            bins.append(bin_min)  # span bins by the the previous bins min + the increment
        return bins

    def avgTest(self, bins, gene, rawdf):
        geneexp = []
        for index, timepoint in enumerate(bins):  # for the index (which is the minimum time) of the bin
            cells_in_this_bin = []
            if index >= len(bins) - 1:  # don't go over range
                pass
            else:
                for i, each in enumerate(rawdf[self.bin_by_var]):
                    # if the cell's ptime is >= the bin, but < bin + 1
                    # then it is within the range/increment and is part of that bin
                    if float(each) >= timepoint and float(each) < bins[index + 1]:
                        cells_in_this_bin.append((i, each))
                # pass the bin and cells in the bin for averaging
                if len(cells_in_this_bin) > 0:
                    # print('bins', timepoint, bins[index + 1])
                    # print('cells in bin ', cells_in_this_bin)
                    geneexp.append(self.lookupCells2(cells_in_this_bin, gene, rawdf))
                else:
                    geneexp.append(np.nan)
        return geneexp

    def maxbinsearch(self):
        maxbin = 5
        bin_range = self.var_max - self.var_min
        weight_dict = {}
        for bin_num in range(4, 101):
            bin_increment = bin_range / bin_num
            bin_min = self.var_min
            bins = [bin_min]
            for i in range(bin_num):
                bin_min += bin_increment
                bins.append(bin_min)
            weights = []
            for index, timepoint in enumerate(bins):  # bin_min, inc1, inc2, inc3, inc4, inc5
                if index >= len(bins) - 1:
                    pass
                else:
                    cells_in_this_bin = []
                    for i, each in enumerate(self.var):
                        # if the cell's ptime is >= the bin, but < bin + 1
                        # then it is within the range/increment and is part of that bin
                        if float(each) >= timepoint and float(each) < bins[index + 1]:
                            cells_in_this_bin.append((i, each))
                    weights.append(len(cells_in_this_bin))
            weight_dict[bin_num] = weights
            for each in weights:
                if each < self.samp_per_bin:
                    if maxbin < 5:
                        return 5, weight_dict
                    else:
                        return maxbin, weight_dict
                    # print('len', bin_num, len(cells_in_this_bin))
                    # if len(cells_in_this_bin) < self.samp_per_bin:
                    #    return maxbin
            maxbin = bin_num
        return maxbin, weight_dict

    def makeSubExp(self, fraction, rawdf, gene, binslist):
        # retrieve a subsample of 80% of the total data
        dfi = rawdf.index.tolist()
        subs = rawdf.iloc[np.random.choice(dfi, math.ceil(len(dfi) * fraction), replace=False)]
        subexp = []
        for bin_i in range(len(binslist) - 1):
            subi = subs.loc[(subs[self.bin_by_var] >= binslist[bin_i]) & (subs[self.bin_by_var] <= binslist[bin_i + 1])]
            if len(subi[gene].tolist()) == 0:
                subexp.append(subexp[-1])
            else:
                subexp.append(np.mean(subi[gene].tolist()))
        return subexp

    def splineParamSelection(self):
        random.seed(42)
        rawexp = pd.DataFrame(self.data_dict[self.branch])
        var_dict = {}
        var_list = []
        for col in rawexp.columns.tolist()[9:]:
            rawexp[col] = pd.to_numeric(rawexp[col])
            var_dict[col] = [rawexp[col].var(), rawexp[col].mean(), rawexp[col].max()]
            var_list.append(rawexp[col].var())
        var_list.sort(reverse=True)
        top100 = var_list[:100]
        varmin = top100[-1]
        genelist = []
        for k, v in var_dict.items():
            if v[0] >= varmin:
                genelist.append(k)
        ############################
        self.var_min = min(self.var)
        self.var_max = max(self.var)
        cellnum = len(self.var)
        if cellnum < (self.samp_per_bin * 4):  # 30 cells per 4 bin minimum
            print('Escape case')
        maxbin, self.weights = self.maxbinsearch()
        smoothing_factors = [0.01, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 7.5, 10]
        sf_by_bin = pd.DataFrame(np.nan, columns=range(4, maxbin), index=smoothing_factors)
        allcoeff = pd.DataFrame(columns=range(4, maxbin), index=smoothing_factors)
        genedict = {}
        for gene in genelist:
            genedf = pd.DataFrame(columns=range(4, maxbin), index=smoothing_factors)  # full of nans
            dfcoeff = pd.DataFrame(columns=range(4, maxbin), index=smoothing_factors)
            genedict[gene] = {}
            for b in range(4, maxbin):  # 101
                genedict[gene][b] = {}
                bins = self.makeBinList(b)
                averaged_exp = self.avgTest(bins, gene, rawexp)
                sf_aic = []
                sf_coeff = []
                for s in smoothing_factors:
                    genedict[gene][b][s] = [], [], [], []
                    spl = interpolate.UnivariateSpline(range(b), averaged_exp, s=s, w=self.weights[
                        b])  # doesn't seem like I can use ptime because of duplicate 0's
                    splineexp = spl(range(b))
                    splineexp[splineexp < 0] = 0
                    knots = len(spl.get_knots())
                    nLL = spl.get_residual()  # should be weighted residuals
                    spl_derivative = spl.derivative()
                    derv = spl_derivative(range(b))
                    rough = derv[-1] - derv[0]  # end - start y derivative values
                    # make sure that ALL genes have an AIC don't want unfair smoothing
                    AIC = (2 * knots) + (2 * nLL) + (2 * s * rough)
                    sf_aic.append(AIC)
                    # coefficient of variation
                    # https://en.wikipedia.org/wiki/Coefficient_of_variation
                    # make subsample splines to test robustness
                    for i in range(30):
                        subexp = self.makeSubExp(0.8, rawexp, gene, bins)
                        subspl = interpolate.UnivariateSpline(range(b), subexp, s=s)
                        subsplineexp = subspl(range(b))
                        subsplineexp[subsplineexp < 0] = 0
                        subcoeffs = subspl.get_coeffs()
                        genedict[gene][b][s][0].append(subcoeffs.tolist())
                        subknots = subspl.get_knots()
                        genedict[gene][b][s][1].append(subknots.tolist())
                        subres = subspl.get_residual()
                        genedict[gene][b][s][2].append(subres)
                        genedict[gene][b][s][3].append(subsplineexp.tolist())
                    d = {}
                    pcoeffvar = []
                    for elem in genedict[gene][b][s][0]:
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
                dfcoeff[b] = sf_coeff
                genedf[b] = sf_aic
            normalized_genedf = (genedf - genedf.min().min()) / (genedf.max().max() - genedf.min().min())
            sf_by_bin = sf_by_bin.add(normalized_genedf, fill_value=0)
            allcoeff = allcoeff.add(dfcoeff, fill_value=0)
        distancedf = pd.DataFrame(columns=sf_by_bin.columns, index=sf_by_bin.index)
        distancedf2 = pd.DataFrame(columns=sf_by_bin.columns, index=sf_by_bin.index)
        writeSpline(sf_by_bin, genedict, allcoeff, self.outdir, self.branch)
        sys.stdout.flush()
        avgaic = sf_by_bin.div(len(genelist))
        avgcoeff = allcoeff.div(len(genelist))
        maxaic = avgaic.max().max() #sf_by_bin.max().max()
        maxcoeff = avgcoeff.max().max() #allcoeff.max().max()
        coeff_thresh = 1
        for col in avgaic.columns:
            for ind in avgaic.index:
                # normalize by the maximum point in the dataframe and calculate the distance from the origin
                # normalization allows for each axis to have an equal effect and distance is used to balance the prefential treatment of one component over another
                distancedf[col][ind] = np.sqrt(((avgaic[col][ind] / maxaic) ** 2) + ((avgcoeff[col][ind] / maxcoeff) ** 2))  # np.sqrt( ((x2 - 0)**2) + ((y2 - 0)**2) )
                if avgcoeff[col][ind] <= coeff_thresh:
                    distancedf2[col][ind] = np.sqrt(((avgaic[col][ind] / maxaic) ** 2) + ((avgcoeff[col][ind] / maxcoeff) ** 2))  # np.sqrt( ((x2 - 0)**2) + ((y2 - 0)**2) )
        print('avgcoeff', self.branch, avgcoeff)
        print('d2', self.branch, distancedf2)
        sys.stdout.flush()
        opt_param = distancedf2.min().min()
        for col in distancedf2.columns.tolist():
            if opt_param in distancedf2[col].tolist() and opt_param >= 0:
                row = distancedf2[col].tolist().index(opt_param)
                print('Branch', self.branch, 'Smoothing Factor of ', smoothing_factors[row], 'Bin number of ', col)
                return smoothing_factors[row], col
        print('Did not pass if conditional!!!')

    def prune_outlier_cells(self):
        exp = pd.DataFrame(self.data_dict[self.branch])
        exp[self.bin_by_var] = [float(i) for i in exp[self.bin_by_var].tolist()]
        Q1, Q3 = np.percentile(exp[self.bin_by_var].tolist(), [25, 75])
        IQR = Q3 - Q1
        ul = Q3 + 1.5 * IQR
        ll = Q1 - 1.5 * IQR
        # upper = exp[exp[self.bin_by_var] > ul]
        exp = exp[(exp[self.bin_by_var] < ul) & (exp[self.bin_by_var] > ll)]
        self.data_dict[self.branch] = exp.to_dict('list')
        return

    def smoother(self):
        # get random set equal or double to largest target set in tf_dict
        # make a factor called 'random' to compare to
        # adjust future calls of random background to reflect this
        self.s_param = None
        self.smoothGeneBinDict = {}
        # write a function to prune outlier cells by ptime
        self.prune_outlier_cells()
        self.var = self.data_dict[self.branch][self.bin_by_var]
        self.samp_per_bin = 10
        if len(self.var) < (self.samp_per_bin * 4):
            print('Not enough cells in this branch to run. Branch = ', self.branch, ' Cellnum = ', len(self.var))
            print('stopping', self.branch)
            return self.smoothGeneBinDict, self.number_of_bins, self.tf_dict, self.s_param
        self.s_param, self.bin_param = self.splineParamSelection()

        self.number_of_bins = self.bin_param
        print('bins after', self.number_of_bins, self.bin_param, self.s_param)

        d_lengths = [len(targets) for targets in self.tf_dict.values()]
        mean_factors = (np.max(d_lengths))
        if mean_factors > len(list(self.data_dict[self.branch].keys())[9:]):
            print("Check the sizing of random and the real data.")
        random.seed(11)
        rTargets1 = random.choices(list(self.data_dict[self.branch].keys())[9:], k=int(np.max(d_lengths)))
        random.seed(22)
        rTargets2 = random.choices(list(self.data_dict[self.branch].keys())[9:], k=int(np.max(d_lengths)))
        random.seed(33)
        rTargets3 = random.choices(list(self.data_dict[self.branch].keys())[9:], k=int(np.max(d_lengths)))
        background1 = [[x, '?'] for x in rTargets1]
        background2 = [[x, '?'] for x in rTargets2]
        background3 = [[x, '?'] for x in rTargets3]
        self.tf_dict['Random1'] = background1
        self.tf_dict['Random2'] = background2
        self.tf_dict['Random3'] = background3
        #### don't exclude Random!!! ####
        self.binningData()  # binning -> getData

        for each_gene in list(self.data_dict[self.branch].keys())[9:]:  # should be only genes
            x = []
            gene_val = []
            for i in range(self.number_of_bins):
                if i in self.avgGeneBinDict:
                    if each_gene in self.avgGeneBinDict[i]:
                        x.append(i)
                        gene_val.append(self.avgGeneBinDict[i][each_gene])
            x_arr = np.array(x)
            g = np.array(gene_val)
            if len(x) >= 4 and len(gene_val) >= 4:
                spl = interpolate.UnivariateSpline(x, gene_val, s=self.s_param, w=self.weights[self.number_of_bins])
                gene_spline = spl(x)
                gene_spline[gene_spline < 0] = 0
                # self.smoothGeneBinDict[each_gene] = gene_spline.tolist()
                minval = min(self.data_dict[self.branch][each_gene])
                m = round(float(minval), 4)
                self.smoothGeneBinDict[each_gene] = [x if x > m else m for x in gene_spline.tolist()]
        for rgene in ['Random1', 'Random2', 'Random3']:
            x = []
            gene_val = []
            for i in range(self.number_of_bins):
                if i in self.avgGeneBinDict:
                    if rgene in self.avgGeneBinDict[i]:
                        x.append(i)
                        gene_val.append(self.avgGeneBinDict[i][rgene])
            x_arr = np.array(x)
            g = np.array(gene_val)
            if len(x) >= 4 and len(gene_val) >= 4:
                spl = interpolate.UnivariateSpline(x, gene_val, s=self.s_param, w=self.weights[self.number_of_bins])
                gene_spline = spl(x_arr)
                gene_spline[gene_spline < 0] = 0
                self.smoothGeneBinDict[rgene] = gene_spline.tolist()

        ptime_med = self.avgGeneBinDict['psuedotime_med']
        print('ptime_med', ptime_med)
        self.smoothGeneBinDict['psuedotime_mean'] = ptime_med
        self.smoothGeneBinDict['psuedotime_med'] = ptime_med
        writeData(self.smoothGeneBinDict, self.branch, 'Smoothed_Expression', '', self.outdir, self.number_of_bins,
                  self.s_param, len(self.var))
        return self.smoothGeneBinDict, self.number_of_bins, self.tf_dict, self.s_param


def writeSpline(aicdf, subsampledict, coeffdf, outdir, branch):
    aicfile = outdir + '/' + '_'.join(branch.split(' -> ')) + '/' + 'aic_df.csv'
    subsampfile = outdir + '/' + '_'.join(branch.split(' -> ')) + '/' + 'subsample_dictionary.txt'
    coefffile = outdir + '/' + '_'.join(branch.split(' -> ')) + '/' + 'coeff_df.csv'
    aicdf.to_csv(aicfile)
    coeffdf.to_csv(coefffile)
    first_line = 'Branch: ' + branch
    subsampf = open(subsampfile, 'x')
    subsampf.write(first_line)
    subsampf.write('\n')
    subsampf.write(str(subsampledict))
    subsampf.close()


def writeData(data, branch, method_type, label, outdir, bins, s_param, cellnum):
    '''
    Written by: Maulding
    '''
    file = outdir + '/' + '_'.join(branch.split(' -> ')) + '/' + label + method_type + '_dictionary.txt'
    f = open(file, 'x')
    first_line = 'Branch: ' + branch + ' Bins: ' + str(bins) + ' Smoothing Factor: ' + str(
        s_param) + ' Cell Count: ' + str(cellnum)
    f.write(first_line)
    f.write('\n')
    f.write(method_type)
    f.write('\n')
    f.write('The following is a dictionary that can be further assessed by the user.')
    f.write('\n')
    f.write(str(data))  # changed to f.write(str(data))
    f.close()


class correlationCalculations:
    '''
    Written by: Maulding
    '''

    def __init__(self, branch, outdir, number_of_bins, bin_by_var, even_distribution, tf_dict, data_dict):
        self.branch = branch
        self.outdir = outdir
        self.number_of_bins = number_of_bins
        self.bin_by_var = bin_by_var
        self.even_distribution = even_distribution
        self.tf_dict = tf_dict
        self.data_dict = data_dict

    def calc_MI(self, x, y, bins):
        c_xy = np.histogram2d(x, y, bins)[0]
        mi = metrics.mutual_info_score(None, None, contingency=c_xy)
        return mi

    def metricRoll(self, fact_name, tar_name, factarr, tararr):
        # select the peak metric based on significance (R, S) and strength (MI, DTW)
        Rcurmax = ((0, 1), 0)  # R, pval, step
        Scurmax = ((0, 1), 0)  # S, pval, step
        MIcurmax = (0, 0)  # MI content, step
        DTWcurmax = (100, 0)  # DTW, step
        R_position_max = []
        S_position_max = []
        MI_position_max = []
        DTW_position_max = []
        for step in range(math.ceil(len(self.smoothGeneBinDict[fact_name]) / 2), len(
                self.smoothGeneBinDict[fact_name])):  # range(3, len(self.smoothGeneBinDict[fact_name]) - 3):
            if np.count_nonzero(factarr[:step]) != 0 and np.count_nonzero(
                    tararr[len(self.smoothGeneBinDict[tar_name]) - step:]) != 0:
                thismax = stats.pearsonr(factarr[:step], tararr[len(self.smoothGeneBinDict[tar_name]) - step:])
                R_position_max.append((thismax, len(self.smoothGeneBinDict[fact_name]) - step))
                if thismax[1] < Rcurmax[0][1]:
                    Rcurmax = (thismax, len(self.smoothGeneBinDict[fact_name]) - step)
                thismax = stats.spearmanr(factarr[:step], tararr[len(self.smoothGeneBinDict[tar_name]) - step:])
                S_position_max.append((thismax, len(self.smoothGeneBinDict[fact_name]) - step))
                if thismax[1] < Scurmax[0][1]:
                    Scurmax = (thismax, len(self.smoothGeneBinDict[fact_name]) - step)
                # mi between fact and tar should be maximized in rollback
                thismax = self.calc_MI(factarr[:step], tararr[len(self.smoothGeneBinDict[tar_name]) - step:],
                                       self.number_of_bins)
                MI_position_max.append((thismax, len(self.smoothGeneBinDict[fact_name]) - step))
                if thismax > MIcurmax[0]:
                    MIcurmax = (thismax, len(self.smoothGeneBinDict[fact_name]) - step)
                # dtw should be minimized
                thismax, cost_matrix, acc_cost_matrix, path = accelerated_dtw(factarr[:step], tararr[len(
                    self.smoothGeneBinDict[tar_name]) - step:],
                                                                              dist='euclidean')  # returns d, cost_matrix, acc_cost_matrix, path
                DTW_position_max.append((thismax, len(self.smoothGeneBinDict[fact_name]) - step))
                if thismax < DTWcurmax[0]:
                    DTWcurmax = (thismax, len(self.smoothGeneBinDict[fact_name]) - step)

        self.rollingfactor2target_Rcors[fact_name][
            tar_name] = Rcurmax  # curmax is the R, pval, and stepback of the most significant correlation
        self.rollingfactor2target_Scors[fact_name][tar_name] = Scurmax
        self.rollingfactor2target_MI[fact_name][tar_name] = MIcurmax
        self.rollingfactor2target_DTW[fact_name][tar_name] = DTWcurmax
        return R_position_max, S_position_max, MI_position_max, DTW_position_max

    def rollCorFigs(self, position_max, factor, target_i, var_min, var_max):
        y1 = []  # R value
        y2 = []  # pvalue
        x = []  # step
        np.seterr(divide='ignore')  # this line turns off possible RuntimeWarnings for y2.append(-np.log10(i[0][1]))
        for i in position_max:
            y1.append(i[0][0] ** 2)
            if i[0][1] == 0:
                y2.append(0)
            else:
                y2.append(-np.log10(i[0][1]))  # RuntimeWarning: divide by zero encountered in log10
            x.append(i[1] - len(self.smoothGeneBinDict[factor]))
        x, y1, y2 = zip(*sorted(zip(x, y1, y2)))
        return [target_i, x, y1, y2, len(self.smoothGeneBinDict[factor]), var_min, var_max]

    def rollOtherFigs(self, position_max, factor, target_i, var_min, var_max):
        y = []  # R value
        x = []  # step
        np.seterr(divide='ignore')  # this line turns off possible RuntimeWarnings for y2.append(-np.log10(i[0][1]))
        for i in position_max:
            y.append(i[0])
            x.append(i[1] - len(self.smoothGeneBinDict[factor]))
        x, y = zip(*sorted(zip(x, y)))
        return [target_i, x, y, len(self.smoothGeneBinDict[factor]), var_min, var_max]

    def roller(self):
        self.smoothGeneBinDict, self.number_of_bins, self.tf_dict, self.s_param = binData(self.branch, self.outdir,
                                                                                          self.number_of_bins,
                                                                                          self.bin_by_var,
                                                                                          self.even_distribution,
                                                                                          self.tf_dict,
                                                                                          self.data_dict).smoother()
        if len(self.smoothGeneBinDict.keys()) == 0:
            print('Nothing in the branch ', self.branch)
            return
        self.rollingfactor2target_Rcors = {}
        self.rollingfactor2target_Scors = {}
        self.rollingfactor2target_MI = {}
        self.rollingfactor2target_DTW = {}
        self.mutualinfo = {}
        self.factor2target_cors = {}
        self.factor2target_scors = {}
        self.factor2target_dtw = {}
        self.regulation = {}
        pearsonFig = {}
        rollingRFig = {}
        rollingSFig = {}
        rollingMIFig = {}
        MIFig = {}
        rollingDTWFig = {}
        dtwFig = {}
        self.Pprune = {}
        self.Sprune = {}
        self.RollPprune = {}
        self.RollSprune = {}
        self.DTWprune = {}
        self.MIprune = {}
        self.RollDTWprune = {}
        self.RollMIprune = {}

        # bin increment for rolling figure generation
        var = self.data_dict[self.branch][self.bin_by_var]  # the vector of the variable as a list
        var_max = max(var)
        var_min = min(var)
        for factor, targets in self.tf_dict.items():
            if factor in self.smoothGeneBinDict:
                self.rollingfactor2target_Rcors[factor] = {}
                self.rollingfactor2target_Scors[factor] = {}
                self.rollingfactor2target_MI[factor] = {}
                self.rollingfactor2target_DTW[factor] = {}
                self.mutualinfo[factor] = {}
                self.factor2target_cors[factor] = {}
                self.factor2target_scors[factor] = {}
                self.factor2target_dtw[factor] = {}
                self.regulation[factor] = {}
                self.Pprune[factor] = {}
                self.Sprune[factor] = {}
                self.RollPprune[factor] = {}
                self.RollSprune[factor] = {}
                self.DTWprune[factor] = {}
                self.MIprune[factor] = {}
                self.RollDTWprune[factor] = {}
                self.RollMIprune[factor] = {}
                rollingRFig[factor] = []
                rollingSFig[factor] = []
                rollingMIFig[factor] = []
                MIFig[factor] = []
                rollingDTWFig[factor] = []
                dtwFig[factor] = []
                heat_arrs = []
                heat_y = []
                roll_arrs = []
                roll_y = []

                nFactor = np.array(self.smoothGeneBinDict[factor])
                # normalize factor -> heat_arrs
                f_heat = nFactor
                fsum = np.nansum(nFactor)
                np.seterr(divide='ignore')  # this line turns off possible RuntimeWarnings
                f_heat2 = [x / fsum if fsum != 0 else 0 for x in
                           f_heat]  # RuntimeWarning: invalid value encountered in double_scalars
                np.seterr(divide='warn')  # this line turns warnings back  on
                heat_arrs.append(f_heat2)
                heat_y.append(factor)

                roll_arrs.append(f_heat2)
                roll_y.append(factor)

                for target_i in targets:
                    if target_i[0] in self.smoothGeneBinDict:
                        ######## data transform ###########

                        nFactor = np.array(self.smoothGeneBinDict[factor])
                        nTarget = np.array(self.smoothGeneBinDict[target_i[0]])

                        if np.count_nonzero(nTarget) == 0 or np.count_nonzero(nFactor) == 0:
                            continue

                        ####### pearson heatmap #########
                        tar_heat = nTarget
                        tarsum = np.nansum(nTarget)
                        # normalize targets tar_heat2
                        tar_heat2 = [x / tarsum for x in tar_heat]
                        heat_arrs.append(tar_heat2)
                        heat_y.append(target_i[0])

                        ###### pearson correlation ########

                        self.factor2target_cors[factor][target_i[0]] = stats.pearsonr(self.smoothGeneBinDict[factor],
                                                                                      self.smoothGeneBinDict[
                                                                                          target_i[0]])

                        ##### spearman correlation ########

                        self.factor2target_scors[factor][target_i[0]] = stats.spearmanr(self.smoothGeneBinDict[factor],
                                                                                        self.smoothGeneBinDict[
                                                                                            target_i[0]])

                        ####### mutual information ########

                        # self.mutualinfo[factor][target_i[0]] = metrics.mutual_info_score(f_heat2, tar_heat2)
                        self.mutualinfo[factor][target_i[0]] = self.calc_MI(f_heat2, tar_heat2,
                                                                            self.number_of_bins)  # updated MI kbased on quantization with 20 bins
                        MIFig[factor].append([target_i, self.mutualinfo[factor][target_i[0]]])

                        ####### rolling correlation #######

                        R_position_max, S_position_max, MI_position_max, DTW_position_max = self.metricRoll(factor,
                                                                                                            target_i[0],
                                                                                                            nFactor,
                                                                                                            nTarget)

                        #### rolling R^2 distribution #####
                        # currently have R_position_max set to make a list of all slices forward and backward
                        if len(R_position_max) != 0:
                            tarFigSpecsR = self.rollCorFigs(R_position_max, factor, target_i, var_min, var_max)
                            rollingRFig[factor].append(tarFigSpecsR)
                        if len(S_position_max) != 0:
                            tarFigSpecsS = self.rollCorFigs(S_position_max, factor, target_i, var_min, var_max)
                            rollingSFig[factor].append(tarFigSpecsS)
                        if len(MI_position_max) != 0:
                            tarFigSpecsMI = self.rollOtherFigs(MI_position_max, factor, target_i, var_min, var_max)
                            rollingMIFig[factor].append(tarFigSpecsMI)
                        if len(DTW_position_max) != 0:
                            tarFigSpecsDTW = self.rollOtherFigs(DTW_position_max, factor, target_i, var_min, var_max)
                            rollingDTWFig[factor].append(tarFigSpecsDTW)
                        np.seterr(divide='warn')  # this line turns warnings back  on

                        ###### dynamic time warping ########
                        norm_targets = np.array(tar_heat2)
                        norm_factor = np.array(f_heat2)
                        d, cost_matrix, acc_cost_matrix, path = accelerated_dtw(norm_targets, norm_factor,
                                                                                dist='euclidean')
                        self.factor2target_dtw[factor][target_i[0]] = d
                        dtwFig[factor].append([target_i, path[0].tolist(), path[1].tolist()])

                        ########### Regulation ##############
                        self.regulation[factor][
                            target_i[0]] = []  # should be Pearson, Spearman, Rolling regulation direction
                        if self.factor2target_cors[factor][target_i[0]][0] < 0:  # if Pearson R value is negative
                            self.regulation[factor][target_i[0]].append('-')
                        else:
                            self.regulation[factor][target_i[0]].append('+')
                        if self.factor2target_scors[factor][target_i[0]][0] < 0:  # if Spearman R value is negative
                            self.regulation[factor][target_i[0]].append('-')
                        else:
                            self.regulation[factor][target_i[0]].append('+')
                        if self.rollingfactor2target_Rcors[factor][target_i[0]][0][
                            0] < 0:  # if Rolling R value is negative
                            self.regulation[factor][target_i[0]].append('-')
                        else:
                            self.regulation[factor][target_i[0]].append('+')
                        if self.rollingfactor2target_Scors[factor][target_i[0]][0][
                            0] < 0:  # if Rolling S value is negative
                            self.regulation[factor][target_i[0]].append('-')
                        else:
                            self.regulation[factor][target_i[0]].append('+')

                        ######## Prune out discordant regulation for non-random entries ########
                        if factor != 'Random':
                            if target_i[1] == self.regulation[factor][target_i[0]][0] or target_i[1] == '?':
                                self.Pprune[factor][target_i[0]] = self.factor2target_cors[factor][target_i[0]]
                            if target_i[1] == self.regulation[factor][target_i[0]][1] or target_i[1] == '?':
                                self.Sprune[factor][target_i[0]] = self.factor2target_scors[factor][target_i[0]]
                            if target_i[1] == self.regulation[factor][target_i[0]][2] or target_i[1] == '?':
                                self.RollPprune[factor][target_i[0]] = self.rollingfactor2target_Rcors[factor][
                                    target_i[0]]
                            if target_i[1] == self.regulation[factor][target_i[0]][3] or target_i[1] == '?':
                                self.RollSprune[factor][target_i[0]] = self.rollingfactor2target_Scors[factor][
                                    target_i[0]]
                            self.DTWprune[factor][target_i[0]] = self.factor2target_dtw[factor][target_i[0]]
                            self.MIprune[factor][target_i[0]] = self.mutualinfo[factor][target_i[0]]
                            self.RollDTWprune[factor][target_i[0]] = self.rollingfactor2target_DTW[factor][target_i[0]]
                            self.RollMIprune[factor][target_i[0]] = self.rollingfactor2target_MI[factor][target_i[0]]

                ##### Pearson / Spearman heatmap #####
                pearsonFig[factor] = [heat_arrs, heat_y, len(self.smoothGeneBinDict[factor])]

        rpearson = []
        rspearman = []
        rrollpearson = []
        rrollspearman = []
        rrollMI = []
        rrollDTW = []
        rDTW = []
        rMI = []
        # print('f2tcors', self.factor2target_cors['Random1'])
        for each_rand_key in ['Random1', 'Random2', 'Random3']:
            rpearson += [(k + '_' + each_rand_key[-1], v) for k, v in self.factor2target_cors[each_rand_key].items()]
            rspearman += [(k + '_' + each_rand_key[-1], v) for k, v in self.factor2target_scors[each_rand_key].items()]
            rrollpearson += [(k + '_' + each_rand_key[-1], v) for k, v in
                             self.rollingfactor2target_Rcors[each_rand_key].items()]
            rrollspearman += [(k + '_' + each_rand_key[-1], v) for k, v in
                              self.rollingfactor2target_Scors[each_rand_key].items()]
            rrollMI += [(k + '_' + each_rand_key[-1], v) for k, v in
                        self.rollingfactor2target_MI[each_rand_key].items()]
            rrollDTW += [(k + '_' + each_rand_key[-1], v) for k, v in
                         self.rollingfactor2target_DTW[each_rand_key].items()]
            rDTW += [(k + '_' + each_rand_key[-1], v) for k, v in self.factor2target_dtw[each_rand_key].items()]
            rMI += [(k + '_' + each_rand_key[-1], v) for k, v in self.mutualinfo[each_rand_key].items()]
            self.factor2target_cors.pop(each_rand_key)
            self.factor2target_scors.pop(each_rand_key)
            self.rollingfactor2target_Rcors.pop(each_rand_key)
            self.rollingfactor2target_Scors.pop(each_rand_key)
            self.rollingfactor2target_MI.pop(each_rand_key)
            self.rollingfactor2target_DTW.pop(each_rand_key)
            self.factor2target_dtw.pop(each_rand_key)
            self.mutualinfo.pop(each_rand_key)
        self.factor2target_cors['Random'] = {i[0]: i[1] for i in rpearson}
        self.factor2target_scors['Random'] = {i[0]: i[1] for i in rspearman}
        self.rollingfactor2target_Rcors['Random'] = {i[0]: i[1] for i in rrollpearson}
        self.rollingfactor2target_Scors['Random'] = {i[0]: i[1] for i in rrollspearman}
        self.rollingfactor2target_MI['Random'] = {i[0]: i[1] for i in rrollMI}
        self.rollingfactor2target_DTW['Random'] = {i[0]: i[1] for i in rrollDTW}
        self.factor2target_dtw['Random'] = {i[0]: i[1] for i in rDTW}
        self.mutualinfo['Random'] = {i[0]: i[1] for i in rMI}

        # Pearson / Spearman heatmap data write
        writeData(pearsonFig, self.branch, 'Correlation', '', self.outdir, self.number_of_bins, self.s_param,
                  len(self.data_dict[self.branch][self.bin_by_var]))
        writeData(rollingRFig, self.branch, 'Rolling_Pearson', '', self.outdir, self.number_of_bins, self.s_param,
                  len(self.data_dict[self.branch][self.bin_by_var]))
        writeData(rollingSFig, self.branch, 'Rolling_Spearman', '', self.outdir, self.number_of_bins, self.s_param,
                  len(self.data_dict[self.branch][self.bin_by_var]))
        writeData(rollingMIFig, self.branch, 'Rolling_MI', '', self.outdir, self.number_of_bins, self.s_param,
                  len(self.data_dict[self.branch][self.bin_by_var]))
        writeData(rollingDTWFig, self.branch, 'Rolling_DTW', '', self.outdir, self.number_of_bins, self.s_param,
                  len(self.data_dict[self.branch][self.bin_by_var]))
        writeData(dtwFig, self.branch, 'DTW', '', self.outdir, self.number_of_bins, self.s_param,
                  len(self.data_dict[self.branch][self.bin_by_var]))
        writeData(MIFig, self.branch, 'Mutual_Information', '', self.outdir, self.number_of_bins, self.s_param,
                  len(self.data_dict[self.branch][self.bin_by_var]))

        return self.rollingfactor2target_Rcors, self.rollingfactor2target_Scors, self.rollingfactor2target_MI, \
               self.rollingfactor2target_DTW, self.mutualinfo, self.factor2target_cors, self.factor2target_scors, \
               self.factor2target_dtw, self.regulation, self.Pprune, self.Sprune, self.RollPprune, self.RollSprune, \
               self.DTWprune, self.MIprune, self.RollDTWprune, self.RollMIprune, self.number_of_bins, self.s_param


def writeResults(data, branch, method_type, outdir):
    '''
    Written by: Maulding
    '''
    file = outdir + '/' + '_'.join(branch.split(' -> ')) + '/' + method_type + '.txt'
    f = open(file, 'x')  # creates new file if it already exists the operation fails
    f.write(branch)
    f.write('\n')
    f.write(method_type)
    sep = "\nFactor\tMean\tStdev\tRandom Mean\tRandom Stdev\tKS Statistic Value\tpvalue\tFDR\t# of Targets\tNames of Targets\tMethod\tBranch\n"
    f.write(sep)
    for each in data:
        newline = str(each[0]) + '\t' + str(each[1]) + '\t' + str(each[2]) + '\t' + str(each[3]) + '\t' + str(
            each[4]) + '\t' + str(each[5]) + '\t' + str(each[6]) + '\t' + str(each[7]) + '\t' + str(
            each[8]) + '\t' + str(each[9]) + '\t' + str(each[10]) + '\t' + branch + '\n'  # changed each[i] to str
        f.write(newline)
    f.close()


class ksStat:
    '''
    Written by: Maulding
    '''

    def __init__(self, branch, outdir, number_of_bins, bin_by_var, even_distribution, tf_dict, data_dict):
        self.branch = branch
        self.outdir = outdir
        self.number_of_bins = number_of_bins
        self.bin_by_var = bin_by_var
        self.even_distribution = even_distribution
        self.tf_dict = tf_dict
        self.data_dict = data_dict
        self.final()

    def loadDictionaries(self):
        try:
            self.rollingfactor2target_Rcors, self.rollingfactor2target_Scors, self.rollingfactor2target_MI, \
            self.rollingfactor2target_DTW, self.mutualinfo, self.factor2target_cors, self.factor2target_scors, \
            self.factor2target_dtw, self.regulation, self.Pprune, self.Sprune, self.RollPprune, self.RollSprune, \
            self.DTWprune, self.MIprune, self.RollDTWprune, self.RollMIprune, self.number_of_bins, self.s_param = \
                correlationCalculations(self.branch, self.outdir, self.number_of_bins, self.bin_by_var,
                                        self.even_distribution, self.tf_dict, self.data_dict).roller()
            return 'go'
        except TypeError:
            print('stop')
            return 'stop'

    def makeCorDistributions(self, non_random, prune):
        self.distribution = {}
        self.prune_distribution = {}
        self.prune_vals = {}
        self.names = {}
        for factor, targets in non_random.items():
            if factor == 'Random':
                rsq_val = []
                for tar, v in targets.items():
                    rsq_val.append((v[0]) ** 2)
                self.distribution[factor] = rsq_val
            else:
                rsquared_vals = []
                rsquared_vals_pruned = []
                self.prune_vals[factor] = {}
                self.names[factor] = []
                for target, value in targets.items():
                    self.names[factor].append(target)
                    rsquared_vals.append((value[0]) ** 2)
                    if target in prune[factor]:
                        rsquared_vals_pruned.append((value[0]) ** 2)
                        self.prune_vals[factor][target] = (value[0]) ** 2
                self.distribution[factor] = rsquared_vals
                self.prune_distribution[factor] = rsquared_vals_pruned

    def makeRollingCorDistribution(self, non_random, prune):
        self.distribution = {}
        self.prune_distribution = {}
        self.prune_vals = {}
        self.names = {}
        for factor, targets in non_random.items():
            if factor == 'Random':
                rand_peak = []
                for tar, v in targets.items():
                    rand_peak.append((v[1]))  # this should be the peak correlation time bin
                self.distribution[factor] = rand_peak
            else:
                peak = []
                prune_peak = []
                self.prune_vals[factor] = {}
                self.names[factor] = []
                for target, value in targets.items():
                    self.names[factor].append(target)
                    peak.append(
                        (value[1]))  # value[0][0] should be R, value[0][1] should be pval, value[1] should be stepback
                    if target in prune[factor]:
                        prune_peak.append((value[1]))
                        self.prune_vals[factor][target] = (value[1])
                self.distribution[factor] = peak
                self.prune_distribution[factor] = prune_peak

    def makeRollingOtherDistribution(self, non_random, prune):
        self.distribution = {}
        self.prune_distribution = {}
        self.prune_vals = {}
        self.names = {}
        for factor, targets in non_random.items():
            if factor == 'Random':
                rand_peak = []
                for tar, v in targets.items():
                    rand_peak.append((v[1]))  # this should be the peak correlation time bin
                self.distribution[factor] = rand_peak
            else:
                peak = []
                prune_peak = []
                self.prune_vals[factor] = {}
                self.names[factor] = []
                for target, value in targets.items():
                    self.names[factor].append(target)
                    peak.append((value[1]))  # value[0] should be MI or DTW, value[1] should be stepback
                    if target in prune[factor]:
                        prune_peak.append((value[1]))
                        self.prune_vals[factor][target] = (value[1])
                self.distribution[factor] = peak
                self.prune_distribution[factor] = prune_peak

    def makeMutualInfoDistribution(self, non_random, prune):
        self.distribution = {}
        self.prune_distribution = {}
        self.prune_vals = {}
        self.names = {}
        for factor, targets in non_random.items():
            if factor == 'Random':
                mutual_val = []
                for tar, v in targets.items():
                    mutual_val.append(v)
                self.distribution[factor] = mutual_val
            else:
                mut_vals = []
                prune_peak = []
                self.prune_vals[factor] = {}
                self.names[factor] = []
                for target, value in targets.items():
                    self.names[factor].append(target)
                    mut_vals.append(value)
                    if target in prune[factor]:
                        prune_peak.append(value)
                        self.prune_vals[factor][target] = value
                self.distribution[factor] = mut_vals
                self.prune_distribution[factor] = prune_peak

    def makeDTWDistribution(self, non_random, prune):
        self.distribution = {}
        self.prune_distribution = {}
        self.prune_vals = {}
        self.names = {}
        for factor, targets in non_random.items():
            if factor == 'Random':
                dtw_val = []
                for tar, v in targets.items():
                    dtw_val.append(v)
                self.distribution[factor] = dtw_val
            else:
                d_vals = []
                prune_peak = []
                self.prune_vals[factor] = {}
                self.names[factor] = []
                for target, value in targets.items():
                    self.names[factor].append(target)
                    d_vals.append(value)
                    if target in prune[factor]:
                        prune_peak.append(value)
                        self.prune_vals[factor][target] = value
                self.distribution[factor] = d_vals
                self.prune_distribution[factor] = prune_peak

    def compareCor(self, method_type):
        self.ks_dict = {}
        self.prune_ks = {}
        good_tars = {}
        for factor, targets in self.tf_dict.items():
            if factor in self.data_dict[self.branch] and factor in self.distribution.keys():
                if len(self.distribution[factor]) != 0 and factor != 'Random' and len(self.distribution['Random']) != 0:
                    self.ks_dict[factor] = stats.ks_2samp(self.distribution[factor],
                                                          self.distribution['Random']), np.mean(
                        self.distribution[factor]), np.std(self.distribution[factor]), np.mean(
                        self.distribution['Random']), np.std(self.distribution['Random']), len(
                        self.distribution[factor]), self.names[factor]
                    # IQR pruning
                    good = []
                    good_tars[factor] = [factor]
                    if len(self.prune_distribution[factor]) > 0:
                        for tar in targets:
                            if tar[0] in self.prune_vals[factor]:
                                if method_type == "Pearson_Correlation" or method_type == "Spearman_Correlation" or method_type == "Mutual_Information":
                                    perc25 = np.percentile(self.prune_distribution[factor], 25)
                                    if self.prune_vals[factor][tar[0]] >= perc25 and len(
                                            self.prune_distribution[factor]) > 4:  # Pearson/Spearman
                                        good.append(self.prune_vals[factor][tar[0]])
                                        good_tars[factor].append(tar[0])
                                    if len(self.prune_distribution[factor]) <= 4:
                                        good.append(self.prune_vals[factor][tar[0]])
                                        good_tars[factor].append(tar[0])
                                if method_type == "Dynamic_Time_Warping":
                                    perc75 = np.percentile(self.prune_distribution[factor], 75)
                                    if self.prune_vals[factor][tar[0]] <= perc75 and len(
                                            self.prune_distribution[factor]) > 4:
                                        good.append(self.prune_vals[factor][tar[0]])
                                        good_tars[factor].append(tar[0])
                                    if len(self.prune_distribution[factor]) <= 4:
                                        good.append(self.prune_vals[factor][tar[0]])
                                        good_tars[factor].append(tar[0])
                                if 'Rolling_' in method_type:  # rolling MI and DTW should be covered here
                                    perc75 = np.percentile(self.prune_distribution[factor], 75)
                                    if self.prune_vals[factor][tar[0]] <= perc75 and len(
                                            self.prune_distribution[factor]) > 4:  # Rolling
                                        good.append(self.prune_vals[factor][tar[0]])
                                        good_tars[factor].append(tar[0])
                                    if len(self.prune_distribution[factor]) <= 4:
                                        good.append(self.prune_vals[factor][tar[0]])
                                        good_tars[factor].append(tar[0])
                        if len(good) >= 3:
                            self.prune_ks[factor] = stats.ks_2samp(good, self.distribution['Random']), np.mean(
                                good), np.std(good), np.mean(self.distribution['Random']), np.std(
                                self.distribution['Random']), len(good), good_tars[factor]
        writeData(good_tars, self.branch, method_type, 'Pruned_', self.outdir, self.number_of_bins, self.s_param,
                  len(self.data_dict[self.branch][self.bin_by_var]))
        self.ksTest(self.ks_dict, method_type)
        method = 'Pruned_' + method_type
        self.ksTest(self.prune_ks, method)

    def ksTest(self, ksData, method_type):
        kslist = []
        p = []
        kspassed = []
        for factor, stat in ksData.items():
            # if number of targets is less than 3 then don't include in fdr testing and append at the bottom
            if stat[5] < 4:
                # factor [0], metric [1], metricSD [2], rand [3], randSD[4], ks val [5], ks pval [6], fdr NA, zscore NA, numTars [7], method [method]
                kspassed.append((factor, stat[1], stat[2], stat[3], stat[4], stat[0][0], stat[0][1], 'NA', stat[5],
                                 ','.join(stat[6]), method_type))
            else:
                kslist.append(
                    (factor, stat[1], stat[2], stat[3], stat[4], stat[0][0], stat[0][1], stat[5], ','.join(stat[6])))
                p.append(stat[0][1])
        if len(p) > 0:
            reject, fdr, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(p, alpha=0.05,
                                                                                             method='fdr_bh',
                                                                                             is_sorted=False,
                                                                                             returnsorted=False)
            kslist2 = []
            for i in range(len(kslist)):
                # factor [0], metric [1], metricSD [2], rand [3], randSD[4], ks [5], pval [6], fdr [fdr[i]], zscore [7], numTars [8], method [method]
                each = (
                kslist[i][0], kslist[i][1], kslist[i][2], kslist[i][3], kslist[i][4], kslist[i][5], kslist[i][6],
                fdr[i], kslist[i][7], kslist[i][8], method_type)
                kslist2.append(each)
                if fdr[i] < 0.05 and each[2] > 0:  # choose significant factors with nonzero deviation
                    if 'Rolling' in each[10]:
                        self.significant.append(each)
                    else:
                        if 'Correlation' in each[10] or 'Mutual_Information' in each[10]:
                            if each[1] > each[3]:
                                self.significant.append(each)
                        if 'Dynamic_Time_Warping' in each[10]:  # and rolling???
                            if each[1] < each[3]:
                                self.significant.append(each)
            kslist2.sort(key=operator.itemgetter(7))
            for factor_item in kspassed:
                kslist2.append(factor_item)
            writeResults(kslist2, self.branch, method_type, self.outdir)  # should write a txt file for the analysis
        else:
            print("Branch empty:", self.branch, "\tmethod:", method_type, len(p))
            sys.stdout.flush()

    def writeRegulation(self, dictionary, branch, outdir):
        self.accuracy = {}
        for factor, targets in self.tf_dict.items():
            if factor in dictionary:
                P = S = RollP = RollS = count = unknown = 0
                for target_i in targets:
                    if target_i[0] in dictionary[factor]:
                        if target_i[1] == '?':
                            unknown += 1
                            pass
                        else:
                            count += 1
                            if dictionary[factor][target_i[0]][0] == target_i[1]:
                                P += 1
                            if dictionary[factor][target_i[0]][1] == target_i[1]:
                                S += 1
                            if dictionary[factor][target_i[0]][2] == target_i[1]:
                                RollP += 1
                            if dictionary[factor][target_i[0]][3] == target_i[1]:
                                RollS += 1
                if count == 0:
                    self.accuracy[factor] = ['na', 'na', 'na', 'na', unknown, unknown + count]
                else:
                    self.accuracy[factor] = [P / count, S / count, RollP / count, RollS / count, unknown,
                                             unknown + count]
        file = outdir + '/' + '_'.join(branch.split(' -> ')) + '/Regulation_dictionary.txt'
        f = open(file, 'x')
        f.write(branch)
        f.write('\n')
        f.write('Regulation')
        f.write('\n')
        f.write('The following is a dictionary that can be further assessed by the user.')
        f.write('\n')
        f.write(str(self.accuracy))
        f.close()

    def final(self):
        move = self.loadDictionaries()
        if move == 'stop':
            print('Nothing in this branch: ', self.branch)
            #os.rmdir(self.outdir + '/' + '_'.join(self.branch.split(' -> ')))
            # OSError: [Errno 39] Directory not empty: 'test20_ebi10/1_2'
            return
        self.significant = []
        self.writeRegulation(self.regulation, self.branch, self.outdir)
        self.makeRollingCorDistribution(self.rollingfactor2target_Rcors, self.RollPprune)
        print("Begin Branch ", self.branch)
        sys.stdout.flush()
        self.compareCor("Rolling_Pearson_Correlation")
        self.makeRollingCorDistribution(self.rollingfactor2target_Scors, self.RollSprune)
        self.compareCor("Rolling_Spearman_Correlation")
        self.makeCorDistributions(self.factor2target_cors, self.Pprune)
        self.compareCor("Pearson_Correlation")
        self.makeCorDistributions(self.factor2target_scors, self.Sprune)
        self.compareCor("Spearman_Correlation")
        self.makeMutualInfoDistribution(self.mutualinfo, self.MIprune)
        self.compareCor("Mutual_Information")
        self.makeDTWDistribution(self.factor2target_dtw, self.DTWprune)
        self.compareCor("Dynamic_Time_Warping")
        self.makeRollingOtherDistribution(self.rollingfactor2target_MI, self.RollMIprune)
        self.compareCor("Rolling_Mutual_Information")
        self.makeRollingOtherDistribution(self.rollingfactor2target_DTW, self.RollDTWprune)
        self.compareCor("Rolling_DTW")
        self.significant.sort(key=operator.itemgetter(5), reverse=True)
        writeResults(self.significant, self.branch, "Summary_of_Significantly_NonRandom_Findings", self.outdir)
        print("Finished Branch ", self.branch)
        sys.stdout.flush()
        return self.number_of_bins, self.s_param


########################################################

def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.

    Written by: Maulding and Meredith
    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else:
        myCommandLine = CommandLine(
            myCommandLine)  # interpret the list passed from the caller of main as the commandline.

    try:

        start_time = time.time()

        pool = Pool(processes=myCommandLine.args.threads)

        # get single-cell JSON data as dictionary
        print("Loading JSON")
        sys.stdout.flush()
        with open(myCommandLine.args.json_dictionary) as j:
            data = json.load(j)[0]
        data_dict_raw = json.loads(data)
        print("JSON loaded")
        # print(data_dict_raw)
        sys.stdout.flush()

        # get factor -> target data as dictionary
        print("Loading tf_dict")
        sys.stdout.flush()
        f = open(myCommandLine.args.tf_dict, 'r')
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

        f2 = open(myCommandLine.args.tf_dict, 'r')
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
        sys.stdout.flush()

        branches = []
        # https://stackoverflow.com/questions/273192/how-can-i-safely-create-a-nested-directory
        os.mkdir(myCommandLine.args.outdir)  # should make trajectory analysis folder
        path = myCommandLine.args.outdir + '/'
        for branch in data_dict_raw:
            os.mkdir(
                os.path.join(path, '_'.join(branch.split(' -> '))))  # should make folders within traj for each branch
            branches.append(branch)
        print("branches:", branches)
        # result = pool.starmap(ksStat, zip(branches, repeat(myCommandLine.args.number_of_bins), repeat(myCommandLine.args.bin_by_var), repeat(myCommandLine.args.even_distribution), repeat(tf_dict_raw), repeat(data_dict_raw)))
        result = pool.starmap(ksStat, zip(branches, repeat(myCommandLine.args.outdir),
                                          repeat(myCommandLine.args.number_of_bins),
                                          repeat(myCommandLine.args.bin_by_var),
                                          repeat(myCommandLine.args.even_distribution), repeat(tfs),
                                          repeat(data_dict_raw)))

        log = open(path + 'run_log.txt', 'w')
        log.write('json_dictionary: ' + str(myCommandLine.args.json_dictionary) + '\n')
        log.write('tf_dict: ' + str(myCommandLine.args.tf_dict) + '\n')
        log.write('outdir: ' + str(myCommandLine.args.outdir) + '\n')
        log.write('bin_by_var: ' + str(myCommandLine.args.bin_by_var) + '\n')
        log.write('DREAMIT2 v1.1')
        log.close()
        print("--- %s seconds ---" % (time.time() - start_time))
        f = open(path + 'time.txt', 'w')
        f.write("--- %s seconds ---" % (time.time() - start_time))
        f.close()

    except Usage as err:
        print(err.msg)


if __name__ == "__main__":
    main()
