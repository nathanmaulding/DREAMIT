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
            sumGene += float(self.data_dict[self.branch][gene][cell_info[0]])
        return sumGene / len(cells_in_bin)


    def avgExpression(self, bin_index, cells_in_this_bin):
        self.avgGeneBinDict[bin_index] = {}
        for factor, targets in self.tf_dict.items():
            if factor in self.data_dict[self.branch]:
                self.avgGeneBinDict[bin_index][factor] = self.lookupCells(cells_in_this_bin, factor)
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
        self.binningData()
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
                    self.smoothGeneBinDict[target_i[0]] = tar_spline,
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


class diffExp:

    def __init__(self, branch, outdir, number_of_bins, bin_by_var, even_distribution, tf_dict, data_dict):
        self.branch = branch
        self.outdir = outdir
        self.number_of_bins = number_of_bins
        self.bin_by_var = bin_by_var
        self.even_distribution = even_distribution
        self.tf_dict = tf_dict
        self.data_dict = data_dict
        self.loadData()

    def loadData(self):
        self.smoothGeneBinDict = binData(self.branch, self.outdir, self.number_of_bins, self.bin_by_var, self.even_distribution,
                                         self.tf_dict, self.data_dict).smoother()
        diff = {}
        num = int(self.number_of_bins / 2)
        for factor, targets in self.tf_dict.items():
            if factor in self.smoothGeneBinDict:
                diff[factor] = []
                fstart = self.smoothGeneBinDict[factor][:num]
                fend = self.smoothGeneBinDict[factor][num:]
                #print(factor, len(fstart), len(fend))
                f = stats.ttest_ind(fstart, fend)
                diff[factor].append(f)
                fchange = (np.mean(fstart) - np.mean(fend)) / np.mean(fstart)
                diff[factor].append(fchange)
                for target_i in targets:
                    if target_i[0] in self.smoothGeneBinDict:
                        diff[target_i[0]] = []
                        tstart = self.smoothGeneBinDict[target_i[0]][:num]
                        tend = self.smoothGeneBinDict[target_i[0]][num:]
                        #print(target_i[0], len(tstart), len(tend))
                        t = stats.ttest_ind(tstart, tend)
                        diff[target_i[0]].append(t)
                        tchange = (np.mean(tstart) - np.mean(tend)) / np.mean(tstart)
                        diff[target_i[0]].append(tchange)
        #print(diff)
        print('here')
        file = self.outdir + '/' + '_'.join(self.branch.split(' -> ')) + '.txt'
        f = open(file, 'x')
        f.write(self.branch)
        f.write('\n')
        f.write("Method")
        sep = '\nGene\tMean\tstat\tpvalue\n'
        f.write(sep)
        for k, v in diff.items():
            newline = k + '\t' + str(v[1]) + '\t' + str(v[0][0]) + '\t' + str(v[0][1]) + '\n'
            f.write(newline)
        f.close()



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
        # branch, number_of_bins, bin_by_var, even_distribution, tf_dict
        # print(myCommandLine.args)  # print the parsed argument string

        pool = Pool(processes=myCommandLine.args.threads)

        # get single-cell JSON data as dictionary
        print("Loading JSON")
        sys.stdout.flush()
        with open(myCommandLine.args.json_dictionary) as j:
            data = json.load(j)[0]
        data_dict_raw = json.loads(data)
        print("JSON loaded")
        sys.stdout.flush()

        # get factor -> target data as dictionary
        print("Loading tf_dict")
        sys.stdout.flush()
        f = open(myCommandLine.args.tf_dict, 'r')
        tfs = {}
        direction_dict = {'Activation': '+', 'Repression': '-', 'Unknown': '?'}
        for eachline in f:
            each = eachline.split()
            if each[0] in tfs:
                tfs[each[0]].append([each[1], direction_dict[each[2]]])
            else:
                tfs[each[0]] = []
                tfs[each[0]].append([each[1], direction_dict[each[2]]])
        print("tf_dict Loaded")
        sys.stdout.flush()

        branches = []
        # https://stackoverflow.com/questions/273192/how-can-i-safely-create-a-nested-directory
        os.mkdir(myCommandLine.args.outdir)  # should make trajectory analysis folder
        path = myCommandLine.args.outdir + '/'
        for branch in data_dict_raw:
            os.mkdir(os.path.join(path, '_'.join(branch.split(' -> '))))  # should make folders within traj for each branch
            branches.append(branch)
        print("branches:", branches)
        #result = pool.starmap(ksStat, zip(branches, repeat(myCommandLine.args.number_of_bins), repeat(myCommandLine.args.bin_by_var), repeat(myCommandLine.args.even_distribution), repeat(tf_dict_raw), repeat(data_dict_raw)))
        result = pool.starmap(diffExp, zip(branches, repeat(myCommandLine.args.outdir), repeat(myCommandLine.args.number_of_bins),
                                          repeat(myCommandLine.args.bin_by_var),
                                          repeat(myCommandLine.args.even_distribution), repeat(tfs),
                                          repeat(data_dict_raw)))


    except Usage as err:
        print(err.msg)

if __name__ == "__main__":
    main()
