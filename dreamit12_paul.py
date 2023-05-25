########### CommandLine ############################

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
        self.parser.add_argument('-number_of_bins', type=int, default=20,
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
        self.smoother()


    def getData(self):
        self.cells = self.data_dict[self.branch]["cell_id"]  # get the index of each cell
        self.var = self.data_dict[self.branch][self.bin_by_var]  # the vector of the variable as a list
        self.cell_var = [{} for _ in self.cells]  # should be cell: var dictionary
        self.var_cell = [{} for _ in self.var]  # should be a var: cell dictionary
        for index, time in enumerate(self.var):  # iter through the list of vars
            self.var[index] = float(time)  # change vars to float
            self.cell_var[index][self.cells[index]] = time  # dict of cells (key) and time (vals)
            self.var_cell[index][self.var[index]] = self.cells[index]  # dict of time (keys) and cells (vals)
        # print(self.var_cell)
        self.var_max = max(self.var)
        self.var_min = min(self.var)


    def lookupCells(self, cells_in_bin, gene):
        # requires the list of cells in the gin and the gene as arguments
        if len(cells_in_bin) == 0:
            return 0
        sumGene = 0
        for cell_info in cells_in_bin:
            # use the index of the cell in question to get the gene expression
            #print("gene", gene, "cell_info", cell_info, "error", self.data_dict[self.branch][gene][cell_info[0]], "sumGene", sumGene, "sumGene type", type(sumGene))
            #print("error type", type(self.data_dict[self.branch][gene][cell_info[0]]))
            sumGene += float(self.data_dict[self.branch][gene][cell_info[0]])
        return sumGene / len(cells_in_bin)


    def avgExpression(self, bin_index, cells_in_this_bin):
        self.avgGeneBinDict[bin_index] = {}
        for factor, targets in self.tf_dict.items():
            if factor in self.data_dict[self.branch]:
                self.avgGeneBinDict[bin_index][factor] = self.lookupCells(cells_in_this_bin, factor)
            if 'Random' in factor:
                rand_index = int(factor[-1]) - 1
                self.avgGeneBinDict[bin_index][factor] = self.lookupCells(cells_in_this_bin, self.randfactors[rand_index])
            for target_i in targets:
                if target_i[0] in self.data_dict[self.branch]:
                    self.avgGeneBinDict[bin_index][target_i[0]] = self.lookupCells(cells_in_this_bin, target_i[0])
        # should return average expression of the gene expression of the cells in that bin
        # shouldn't even need to return it #return self.avgGeneBinDict


    def binningData(self):
        self.getData()
        self.avgGeneBinDict = {}
        bin_range = self.var_max - self.var_min  # get min and max time on branch
        bin_increment = bin_range / self.number_of_bins  # create an increment to move by
        bin_min = self.var_min
        bins = [self.var_min]  # create bins based on minimum increments
        self.bin_dict = []
        self.randfactors = random.choices([k for k in self.tf_dict.keys() if k in list(self.data_dict[self.branch].keys())[9:]], k=3)
        for i in range(self.number_of_bins):
            self.bin_dict.append({})
            bin_min += bin_increment
            bins.append(bin_min)   # span bins by the the previous bins min + the increment
        for index, timepoint in enumerate(bins):  # for the index (which is the minimum time) of the bin
            cells_in_this_bin = []
            if index >= len(bins) - 1:  # don't go over range
                pass
            else:
                for i, each in enumerate(self.data_dict[self.branch][self.bin_by_var]):
                    # if the cell's ptime is >= the bin, but < bin + 1
                    # then it is within the range/increment and is part of that bin
                    if float(each) >= timepoint and float(each) < bins[index + 1]:
                        cells_in_this_bin.append((i, each))
                # pass the bin and cells in the bin for averaging
                self.avgExpression(index, cells_in_this_bin)
        return self.avgGeneBinDict


    def smoother(self):
        # get random set equal or double to largest target set in tf_dict
        # make a factor called 'random' to compare to
        # adjust future calls of random background to reflect this
        d_lengths = [len(targets) for targets in self.tf_dict.values()]
        mean_factors = (np.max(d_lengths))
        if mean_factors > len(list(self.data_dict[self.branch].keys())[9:]):
            print("Check the sizing of random and the real data.")
        rTargets1 = random.choices(list(self.data_dict[self.branch].keys())[9:], k=int(np.max(d_lengths)))
        rTargets2 = random.choices(list(self.data_dict[self.branch].keys())[9:], k=int(np.max(d_lengths)))
        rTargets3 = random.choices(list(self.data_dict[self.branch].keys())[9:], k=int(np.max(d_lengths)))
        background1 = [[x, '?'] for x in rTargets1]
        background2 = [[x, '?'] for x in rTargets2]
        background3 = [[x, '?'] for x in rTargets3]
        self.tf_dict['Random1'] = background1
        self.tf_dict['Random2'] = background2
        self.tf_dict['Random3'] = background3
        #### don't exclude Random!!! ####
        self.binningData()  # binning -> getData
        self.smoothGeneBinDict = {}
        for factor, targets in self.tf_dict.items():
            for target_i in targets:
                x = []
                t_val = []
                for i in range(self.number_of_bins):
                    if i in self.avgGeneBinDict:
                        if target_i[0] in self.avgGeneBinDict[i]:
                            x.append(i)
                            t_val.append(self.avgGeneBinDict[i][target_i[0]])
                x_arr = np.array(x)
                t = np.array(t_val)
                if len(x) >= 4 and len(t_val) >= 4:
                    spl = interpolate.UnivariateSpline(x, t_val)
                if len(x_arr) >= 4:
                    tar_spline = spl(x_arr)
                    tar_spline[tar_spline < 0] = 0
                    self.smoothGeneBinDict[target_i[0]] = tar_spline
            bins = []
            f_val = []
            for index in range(self.number_of_bins):
                if index in self.avgGeneBinDict:
                    if factor in self.avgGeneBinDict[index]:
                        bins.append(index)
                        f_val.append(self.avgGeneBinDict[index][factor])
            bins_arr = np.array(bins)
            f = np.array(f_val)
            if len(bins) >= 4 and len(f_val) >= 4:
                fspl = interpolate.UnivariateSpline(bins, f_val)
            if len(bins_arr) >= 4:
                factor_spline = fspl(bins_arr)
                factor_spline[factor_spline < 0] = 0
                self.smoothGeneBinDict[factor] = factor_spline
        return self.smoothGeneBinDict

'''
class randomCompare:

    # Written by: Maulding

    def __init__(self, branch, outdir, number_of_bins, bin_by_var, even_distribution, tf_dict, data_dict):
        self.branch = branch
        self.outdir = outdir
        self.number_of_bins = number_of_bins
        self.bin_by_var = bin_by_var
        self.even_distribution = even_distribution
        self.tf_dict = tf_dict
        self.data_dict = data_dict
        self.randomize()

    def randomize(self):
        d_lengths = [len(targets) for targets in self.tf_dict.values()]
        mean_factors = (np.max(d_lengths)) * 2
        if mean_factors > len(list(self.data_dict[self.branch].keys())[9:]):
            assert "mean_factors variable for selecting random background is too large"
        self.randDict = {}
        for factor, targets in self.tf_dict.items():
            rTargets = random.choices(list(self.data_dict[self.branch].keys())[9:], k=int(mean_factors))
            rTargets = [[r, '?'] for r in rTargets]
            self.randDict[factor] = rTargets
        return self.randDict
'''

def writeData(data, branch, method_type, label, outdir):
    '''
    Written by: Maulding
    '''
    file = outdir + '/' + '_'.join(branch.split(' -> ')) + '/' + label + method_type + '_dictionary.txt'
    f = open(file, 'x')
    f.write(branch)
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
        Rcurmax = ((0, 1), 0)
        Scurmax = ((0, 1), 0)
        MIcurmax = (100, 0)
        DTWcurmax = (100, 0)
        R_position_max = []
        S_position_max = []
        MI_position_max = []
        DTW_position_max = []
        for step in range(3, len(self.smoothGeneBinDict[fact_name]) - 3):
            if np.count_nonzero(factarr[:step]) != 0 and np.count_nonzero(tararr[len(self.smoothGeneBinDict[tar_name]) - step:]) != 0:
                thismax = stats.pearsonr(factarr[:step], tararr[len(self.smoothGeneBinDict[tar_name]) - step:])
                R_position_max.append((thismax, step))
                if thismax[1] < Rcurmax[0][1]:
                    Rcurmax = (thismax, step)
                thismax = stats.spearmanr(factarr[:step], tararr[len(self.smoothGeneBinDict[tar_name]) - step:])
                S_position_max.append((thismax, step))
                if thismax[1] < Scurmax[0][1]:
                    Scurmax = (thismax, step)
                thismax = self.calc_MI(factarr[:step], tararr[len(self.smoothGeneBinDict[tar_name]) - step:], self.number_of_bins)
                MI_position_max.append((thismax, step))
                if thismax < MIcurmax[0]:
                    MIcurmax = (thismax, step)
                thismax, cost_matrix, acc_cost_matrix, path = accelerated_dtw(factarr[:step], tararr[len(self.smoothGeneBinDict[tar_name]) - step:], dist='euclidean') # returns d, cost_matrix, acc_cost_matrix, path
                DTW_position_max.append((thismax, step))
                if thismax < DTWcurmax[0]:
                    DTWcurmax = (thismax, step + len(self.smoothGeneBinDict[tar_name]))
            if np.count_nonzero(factarr[len(self.smoothGeneBinDict[fact_name]) - step:]) != 0 and np.count_nonzero(tararr[:step]) != 0:
                thismax = stats.pearsonr(factarr[len(self.smoothGeneBinDict[fact_name]) - step:], tararr[:step])
                R_position_max.append((thismax, -(step)))
                if thismax[1] < Rcurmax[0][1]:
                    Rcurmax = (thismax, step + len(self.smoothGeneBinDict[tar_name]))
                thismax = stats.spearmanr(factarr[len(self.smoothGeneBinDict[fact_name]) - step:], tararr[:step])
                S_position_max.append((thismax, -(step)))
                if thismax[1] < Scurmax[0][1]:
                    Scurmax = (thismax, step + len(self.smoothGeneBinDict[tar_name]))
                thismax = self.calc_MI(factarr[len(self.smoothGeneBinDict[fact_name]) - step:], tararr[:step], self.number_of_bins)
                MI_position_max.append((thismax, -(step)))
                if thismax < MIcurmax[0]:
                    MIcurmax = (thismax, step + len(self.smoothGeneBinDict[tar_name]))
                thismax, cost_matrix, acc_cost_matrix, path = accelerated_dtw(factarr[len(self.smoothGeneBinDict[fact_name]) - step:], tararr[:step], dist='euclidean') # returns d, cost_matrix, acc_cost_matrix, path
                DTW_position_max.append((thismax, -(step)))
                if thismax < DTWcurmax[0]:
                    DTWcurmax = (thismax, step + len(self.smoothGeneBinDict[tar_name]))
        self.rollingfactor2target_Rcors[fact_name][tar_name] = Rcurmax  # curmax is the R, pval, and stepback of the most significant correlation
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
        self.smoothGeneBinDict = binData(self.branch, self.outdir, self.number_of_bins, self.bin_by_var, self.even_distribution, self.tf_dict, self.data_dict).smoother()
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
            if factor in self.smoothGeneBinDict or 'Random' in factor:
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
                f_heat2 = [x / fsum if fsum != 0 else 0 for x in f_heat]  # RuntimeWarning: invalid value encountered in double_scalars
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

                        self.factor2target_cors[factor][target_i[0]] = stats.pearsonr(self.smoothGeneBinDict[factor], self.smoothGeneBinDict[target_i[0]])

                        ##### spearman correlation ########

                        self.factor2target_scors[factor][target_i[0]] = stats.spearmanr(self.smoothGeneBinDict[factor], self.smoothGeneBinDict[target_i[0]])

                        ####### mutual information ########

                        # self.mutualinfo[factor][target_i[0]] = metrics.mutual_info_score(f_heat2, tar_heat2)
                        self.mutualinfo[factor][target_i[0]] = self.calc_MI(f_heat2, tar_heat2, self.number_of_bins) # updated MI kbased on quantization with 20 bins
                        MIFig[factor].append([target_i, self.mutualinfo[factor][target_i[0]]])

                        ####### rolling correlation #######

                        R_position_max, S_position_max, MI_position_max, DTW_position_max = self.metricRoll(factor, target_i[0], nFactor, nTarget)

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
                        np.seterr(divide='warn') # this line turns warnings back  on

                        ###### dynamic time warping ########
                        norm_targets = np.array(tar_heat2)
                        norm_factor = np.array(f_heat2)
                        d, cost_matrix, acc_cost_matrix, path = accelerated_dtw(norm_targets, norm_factor, dist='euclidean')
                        self.factor2target_dtw[factor][target_i[0]] = d
                        dtwFig[factor].append([target_i, path[0].tolist(), path[1].tolist()])

                        ########### Regulation ##############
                        self.regulation[factor][target_i[0]] = []  # should be Pearson, Spearman, Rolling regulation direction
                        if self.factor2target_cors[factor][target_i[0]][0] < 0:  # if Pearson R value is negative
                            self.regulation[factor][target_i[0]].append('-')
                        else:
                            self.regulation[factor][target_i[0]].append('+')
                        if self.factor2target_scors[factor][target_i[0]][0] < 0:  # if Spearman R value is negative
                            self.regulation[factor][target_i[0]].append('-')
                        else:
                            self.regulation[factor][target_i[0]].append('+')
                        if self.rollingfactor2target_Rcors[factor][target_i[0]][0][0] < 0:  # if Rolling R value is negative
                            self.regulation[factor][target_i[0]].append('-')
                        else:
                            self.regulation[factor][target_i[0]].append('+')
                        if self.rollingfactor2target_Scors[factor][target_i[0]][0][0] < 0:  # if Rolling S value is negative
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
                                self.RollPprune[factor][target_i[0]] = self.rollingfactor2target_Rcors[factor][target_i[0]]
                            if target_i[1] == self.regulation[factor][target_i[0]][3] or target_i[1] == '?':
                                self.RollSprune[factor][target_i[0]] = self.rollingfactor2target_Scors[factor][target_i[0]]
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
        for each_rand_key in ['Random1', 'Random2', 'Random3']:
            rpearson += [(k + '_' + each_rand_key[-1], v) for k, v in self.factor2target_cors[each_rand_key].items()]
            rspearman += [(k + '_' + each_rand_key[-1], v) for k, v in self.factor2target_scors[each_rand_key].items()]
            rrollpearson += [(k + '_' + each_rand_key[-1], v) for k, v in self.rollingfactor2target_Rcors[each_rand_key].items()]
            rrollspearman += [(k + '_' + each_rand_key[-1], v) for k, v in self.rollingfactor2target_Scors[each_rand_key].items()]
            rrollMI += [(k + '_' + each_rand_key[-1], v) for k, v in self.rollingfactor2target_MI[each_rand_key].items()]
            rrollDTW += [(k + '_' + each_rand_key[-1], v) for k, v in self.rollingfactor2target_DTW[each_rand_key].items()]
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
        writeData(pearsonFig, self.branch, 'Correlation', '', self.outdir)
        writeData(rollingRFig, self.branch, 'Rolling_Pearson', '', self.outdir)
        writeData(rollingSFig, self.branch, 'Rolling_Spearman', '', self.outdir)
        writeData(rollingMIFig, self.branch, 'Rolling_MI', '', self.outdir)
        writeData(rollingDTWFig, self.branch, 'Rolling_DTW', '', self.outdir)
        writeData(dtwFig, self.branch, 'DTW', '', self.outdir)
        writeData(MIFig, self.branch, 'Mutual_Information', '', self.outdir)

        return self.rollingfactor2target_Rcors, self.rollingfactor2target_Scors, self.rollingfactor2target_MI, \
                   self.rollingfactor2target_DTW, self.mutualinfo, self.factor2target_cors, self.factor2target_scors, \
                   self.factor2target_dtw, self.regulation, self.Pprune, self.Sprune, self.RollPprune, self.RollSprune, \
                   self.DTWprune, self.MIprune, self.RollDTWprune, self.RollMIprune

class smoothGenes:
    '''This function flows from init -> smoother -> binningData[bin] -> getData -> avgExpression
        -> lookupCells -> make avgGeneDict binningData[bin] -> next bin -> smoother

        Written by: Maulding
        '''

    def __init__(self, branch, outdir, number_of_bins, bin_by_var, even_distribution, data_dict):
        self.branch = branch
        self.outdir = outdir
        self.number_of_bins = number_of_bins
        self.bin_by_var = bin_by_var
        self.even_distribution = even_distribution
        self.data_dict = data_dict
        self.smooth()

    def cellData(self):
        self.cells = self.data_dict[self.branch]["cell_id"]  # get the index of each cell
        self.var = self.data_dict[self.branch][self.bin_by_var]  # the vector of the variable as a list
        self.cell_var = [{} for _ in self.cells]  # should be cell: var dictionary
        self.var_cell = [{} for _ in self.var]  # should be a var: cell dictionary
        for index, time in enumerate(self.var):  # iter through the list of vars
            self.var[index] = float(time)  # change vars to float
            self.cell_var[index][self.cells[index]] = time  # dict of cells (key) and time (vals)
            self.var_cell[index][self.var[index]] = self.cells[index]  # dict of time (keys) and cells (vals)
        # print(self.var_cell)
        self.var_max = max(self.var)
        self.var_min = min(self.var)

    def lookup(self, cells_in_bin, gene):
        # requires the list of cells in the gin and the gene as arguments
        if len(cells_in_bin) == 0:
            return 0
        sumGene = 0
        for cell_info in cells_in_bin:
            # use the index of the cell in question to get the gene expression
            #print("gene", gene, "cell_info", cell_info, "error", cell_info[0], "sumGene", sumGene, "sumGene type", type(sumGene))
            # print("error type", type(self.data_dict[self.branch][gene][cell_info[0]]))
            sumGene += float(self.data_dict[self.branch][gene][cell_info[0]])
        return sumGene / len(cells_in_bin)

    def averageExp(self, bin_i, cells_in_bin):
        self.averaged[bin_i] = {}
        for each_gene in list(self.data_dict[self.branch].keys())[9:]: # should be only genes
            #print("should be gene", each_gene)
            self.averaged[bin_i][each_gene] = self.lookup(cells_in_bin, each_gene)
        # should return average expression of the gene expression of the cells in that bin
        # shouldn't even need to return it #return self.avgGeneBinDict

    def bin(self):
        self.cellData()
        self.averaged = {}
        bin_range = self.var_max - self.var_min  # get min and max time on branch
        bin_increment = bin_range / self.number_of_bins  # create an increment to move by
        bin_min = self.var_min
        bins = [self.var_min]  # create bins based on minimum increments
        self.bin_dict = []
        for i in range(self.number_of_bins):
            self.bin_dict.append({})
            bin_min += bin_increment
            bins.append(bin_min)  # span bins by the the previous bins min + the increment
        for index, timepoint in enumerate(bins):  # for the index (which is the minimum time) of the bin
            cells_in_this_bin = []
            if index >= len(bins) - 1:  # don't go over range
                pass
            else:
                for i, each in enumerate(self.data_dict[self.branch][self.bin_by_var]):
                    # if the cell's ptime is >= the bin, but < bin + 1
                    # then it is within the range/increment and is part of that bin
                    if float(each) >= timepoint and float(each) < bins[index + 1]:
                        cells_in_this_bin.append((i, each))
                # pass the bin and cells in the bin for averaging
                self.averageExp(index, cells_in_this_bin)
        return self.averaged

    def smooth(self):
        self.bin()
        self.smoothgenes = {}
        for each_gene in list(self.data_dict[self.branch].keys())[9:]:  # should be only genes
            x = []
            gene_val = []
            for i in range(self.number_of_bins):
                if i in self.averaged:
                    if each_gene in self.averaged[i]:
                        x.append(i)
                        gene_val.append(self.averaged[i][each_gene])
            x_arr = np.array(x)
            g = np.array(gene_val)
            if len(x) >= 4 and len(gene_val) >= 4:
                spl = interpolate.UnivariateSpline(x, gene_val)
                gene_spline = spl(x_arr)
                gene_spline[gene_spline < 0] = 0
                self.smoothgenes[each_gene] = gene_spline.tolist()
        writeData(self.smoothgenes, self.branch, 'Smoothed_Expression', '', self.outdir)


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
        newline = str(each[0]) + '\t' + str(each[1]) + '\t' + str(each[2]) + '\t' + str(each[3]) + '\t' + str(each[4]) + '\t' + str(each[5]) + '\t' + str(each[6]) + '\t' + str(each[7]) + '\t' + str(each[8]) + '\t' + str(each[9]) + '\t' + str(each[10]) + '\t' + branch + '\n'  # changed each[i] to str
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
        self.rollingfactor2target_Rcors, self.rollingfactor2target_Scors, self.rollingfactor2target_MI, \
        self.rollingfactor2target_DTW, self.mutualinfo, self.factor2target_cors, self.factor2target_scors, \
        self.factor2target_dtw, self.regulation, self.Pprune, self.Sprune, self.RollPprune, self.RollSprune, \
        self.DTWprune, self.MIprune, self.RollDTWprune, self.RollMIprune = \
            correlationCalculations(self.branch, self.outdir, self.number_of_bins, self.bin_by_var, self.even_distribution, self.tf_dict, self.data_dict).roller()
        smoothGenes(self.branch, self.outdir, self.number_of_bins, self.bin_by_var, self.even_distribution, self.data_dict)

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
        print('rollP', non_random.keys())
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
                    peak.append((value[1]))  # value[0][0] should be R, value[0][1] should be pval, value[1] should be stepback
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
                if len(self.distribution[factor]) != 0 and factor != 'Random':
                    # eliminated zobs from index 5
                    #print(method_type, self.distribution.keys())
                    #print(self.distribution['Random'])
                    self.ks_dict[factor] = stats.ks_2samp(self.distribution[factor], self.distribution['Random']), np.mean(self.distribution[factor]), np.std(self.distribution[factor]), np.mean(self.distribution['Random']), np.std(self.distribution['Random']), len(self.distribution[factor]), self.names[factor]
                    # IQR pruning
                    good = []
                    good_tars[factor] = [factor]
                    if len(self.prune_distribution[factor]) > 0:
                        # mean_std = [np.mean(self.prune_distribution[factor]), np.std(self.prune_distribution[factor])]
                        for tar in targets:
                            if tar[0] in self.prune_vals[factor]:
                                if method_type == "Pearson_Correlation" or method_type == "Spearman_Correlation" or method_type == "Mutual_Information":
                                    perc25 = np.percentile(self.prune_distribution[factor], 25)
                                    # if self.prune_vals[factor][tar[0]] >= (mean_std[0] - mean_std[1]) and len(self.prune_distribution[factor]) > 4:
                                    if self.prune_vals[factor][tar[0]] >= perc25 and len(self.prune_distribution[factor]) > 4: # Pearson/Spearman
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
                                if 'Rolling_' in method_type: # rolling MI and DTW should be covered here
                                    perc75 = np.percentile(self.prune_distribution[factor], 75)
                                    # if self.prune_vals[factor][tar[0]] <= (mean_std[0] + mean_std[1]) and len(self.prune_distribution[factor]) > 4:  # Rolling
                                    if self.prune_vals[factor][tar[0]] <= perc75 and len(self.prune_distribution[factor]) > 4:  # Rolling
                                        good.append(self.prune_vals[factor][tar[0]])
                                        good_tars[factor].append(tar[0])
                                    if len(self.prune_distribution[factor]) <= 4:
                                        good.append(self.prune_vals[factor][tar[0]])
                                        good_tars[factor].append(tar[0])
                        if len(good) >= 3:
                            self.prune_ks[factor] = stats.ks_2samp(good, self.distribution['Random']), np.mean(good), np.std(good), np.mean(self.distribution['Random']), np.std(self.distribution['Random']), len(good), good_tars[factor]
        writeData(good_tars, self.branch, method_type, 'Pruned_', self.outdir)
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
                kspassed.append((factor, stat[1], stat[2], stat[3], stat[4], stat[0][0], stat[0][1], 'NA', stat[5], ','.join(stat[6]), method_type))
            else:
                kslist.append((factor, stat[1], stat[2], stat[3], stat[4], stat[0][0], stat[0][1], stat[5], ','.join(stat[6])))
                p.append(stat[0][1])
        if len(p) > 0:
            reject, fdr, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(p, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
            kslist2 = []
            for i in range(len(kslist)):
                # factor [0], metric [1], metricSD [2], rand [3], randSD[4], ks [5], pval [6], fdr [fdr[i]], zscore [7], numTars [8], method [method]
                each = (kslist[i][0], kslist[i][1], kslist[i][2], kslist[i][3], kslist[i][4], kslist[i][5], kslist[i][6], fdr[i], kslist[i][7], kslist[i][8], method_type)
                kslist2.append(each)
                if fdr[i] < 0.05:
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
                    self.accuracy[factor] = [P/count, S/count, RollP/count, RollS/count, unknown, unknown + count]
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
        self.loadDictionaries()
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


########################################################

def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.

    Written by: Maulding and Meredith
    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else:
        myCommandLine = CommandLine(myCommandLine)  # interpret the list passed from the caller of main as the commandline.

    try:

        start_time = time.time()

        pool = Pool(processes=myCommandLine.args.threads)

        # get single-cell JSON data as dictionary
        print("Loading JSON")
        sys.stdout.flush()
        with open(myCommandLine.args.json_dictionary) as j:
            data_dict_raw = json.load(j)
        #data_dict_raw = json.loads(data)
        print("JSON loaded")
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
                if each[1] in doubles[each[0]]:
                    if [each[1], '?'] not in tfs[each[0]]:
                        tfs[each[0]].append([each[1], '?'])
                else:
                    tfs[each[0]].append([each[1], direction_dict[each[2]]])
            else:
                tfs[each[0]] = []
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
        for branch in data_dict_raw['branches']:
            os.mkdir(os.path.join(path, '_'.join(branch.split(' -> '))))  # should make folders within traj for each branch
            branches.append(branch)
        print("branches:", branches)
        #result = pool.starmap(ksStat, zip(branches, repeat(myCommandLine.args.number_of_bins), repeat(myCommandLine.args.bin_by_var), repeat(myCommandLine.args.even_distribution), repeat(tf_dict_raw), repeat(data_dict_raw)))
        result = pool.starmap(ksStat, zip(branches, repeat(myCommandLine.args.outdir), repeat(myCommandLine.args.number_of_bins),
                                          repeat(myCommandLine.args.bin_by_var),
                                          repeat(myCommandLine.args.even_distribution), repeat(tfs),
                                          repeat(data_dict_raw['branches'])))

        print("--- %s seconds ---" % (time.time() - start_time))
        f = open(path + 'time.txt', 'w')
        f.write("--- %s seconds ---" % (time.time() - start_time))
        f.close()

    except Usage as err:
        print(err.msg)

if __name__ == "__main__":
    main()