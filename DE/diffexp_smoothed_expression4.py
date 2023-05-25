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
        self.parser.add_argument('-expression_path',
                                 help="path to directory of branches; expression files are contained within the branch directory")
        self.parser.add_argument('-json_dictionary', help="path to raw json dictionary")
        self.parser.add_argument('-threads', type=int, default=1,
                                 help='input the number of threads to use. Default = 1')

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

import sys
import os
from scipy import stats
import numpy as np
import statsmodels.stats.multitest
import operator
from multiprocessing import Pool
import time
import pandas as pd
import json
from itertools import repeat


# folder = '/Users/nathandmaulding/Desktop/DREAMIT2/EBI-4-Induced dendritic/diffexp/'

def openFile(path):
    data = open(path, 'r')
    figDict = ''
    branch = ''
    for i, line in enumerate(data):
        if i == 0:
            b = line.strip('\n')[8:14]
            branch = '_'.join(b.split(' -> '))
        if i >= 3:
            figDict = figDict + line
    figDict = eval(figDict.replace('inf', 'float("inf")'))  # NameError: name 'inf' not defined (EBI6)
    return figDict, branch


def openJson(path, ivar):
    data = open(path, 'r')
    expDict = ''
    for i, line in enumerate(data):
        line.split()
        if i == ivar:
            expDict = line
    expDict = eval(expDict)
    return expDict


def writeResults(data, branch, method_type, outdir):
    '''
    Written by: Maulding
    '''
    # dir/1_2_diffexp.txt
    file = outdir + '_'.join(branch.split(' -> ')) + '_' + method_type + '.txt'
    f = open(file, 'x')  # creates new file if it already exists the operation fails
    f.write(branch)
    f.write('\n')
    f.write(method_type)
    sep = "\nGene\tFold Change\tTstat\tpvalue\tFDR\n"
    f.write(sep)
    for each in data:
        newline = str(each[0]) + '\t' + str(each[1]) + '\t' + str(each[2]) + '\t' + str(each[3]) + '\t' + str(
            each[4]) + '\n'  # changed each[i] to str
        f.write(newline)
    f.close()


class differential_genes:

    def __init__(self, expression_path, empty):
        self.expression_path = expression_path
        self.empty = empty
        self.diffexp_run()

    def diffexp_run(self):

        folder = self.expression_path + '/'

        for file in os.listdir(folder):
            if 'Smoothed_Expression_dictionary' not in file:
                pass
            else:
                path = folder + file
                exp, branch = openFile(path)
                diff_results = []
                ps = []
                for each_gene in exp.keys():
                    if each_gene in ['psuedotime_mean', 'psuedotime_med', 'cell_ids']:
                        continue
                    half = int(len(exp[each_gene]) / 2)
                    start = exp[each_gene][:half]
                    end = exp[each_gene][half:]
                    # print('start end', type(start), type(end))
                    # print('start end', type(start[0]), type(end[0]))
                    s = stats.ttest_ind(start, end)
                    # print('s', s)
                    pval = s[1]
                    teststat = s[0]
                    diff = np.mean(end) - np.mean(start)
                    if pval == pval:
                        diff_results.append([each_gene, diff, teststat, pval])
                        ps.append(pval)
                if len(ps) > 0:
                    reject, fdr, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(ps, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
                    for i in range(len(diff_results)):
                        diff_results[i].append(fdr[i])
                    diff_results.sort(key=operator.itemgetter(4))
                    writeResults(diff_results, branch, 'smooth_diffexp', folder)

class first_and_last:

    def __init__(self, expression_path, empty):
        self.expression_path = expression_path
        self.empty = empty
        self.diffexp_run()

    def diffexp_run(self):

        folder = self.expression_path + '/'

        for file in os.listdir(folder):
            if 'Smoothed_Expression_dictionary' not in file:
                pass
            else:
                path = folder + file
                exp, branch = openFile(path)
                diff_results = []
                ps = []
                for each_gene in exp.keys():
                    if each_gene in ['psuedotime_mean', 'psuedotime_med', 'cell_ids']:
                        continue
                    start = exp[each_gene][0]
                    end = exp[each_gene][-1]
                    # print('start end', type(start), type(end))
                    # print('start end', type(start[0]), type(end[0]))
                    s = stats.ttest_ind(start, end)
                    # print('s', s)
                    pval = s[1]
                    teststat = s[0]
                    diff = np.mean(end) - np.mean(start)
                    if pval == pval:
                        diff_results.append([each_gene, diff, teststat, pval])
                        ps.append(pval)
                if len(ps) > 0:
                    reject, fdr, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(ps, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
                    for i in range(len(diff_results)):
                        diff_results[i].append(fdr[i])
                    diff_results.sort(key=operator.itemgetter(4))
                    writeResults(diff_results, branch, 'firstlast_diffexp', folder)


class rawdiff:

    def __init__(self, expression_path, json_dictionary):
        self.expression_path = expression_path
        self.json_dictionary = json_dictionary
        self.raw_diffexp()

    def raw_diffexp(self):
        with open(self.json_dictionary) as j:
            data = json.load(j)[0]  # omit [0] for paul
        exp = json.loads(data)  # omit for paul
        brs = exp.keys()
        var = 'psuedotime'
        branches = [b for b in brs if '_'.join(b.split(' -> ')) in os.listdir(self.expression_path)]
        for branch in branches:
            rawexp = pd.DataFrame(exp[branch])
            # prune outlier cells
            rawexp[var] = [float(i) for i in rawexp[var].tolist()]
            Q1, Q3 = np.percentile(rawexp[var], [25, 75])
            IQR = Q3 - Q1
            ul = Q3 + 1.5 * IQR
            ll = Q1 - 1.5 * IQR
            rawexp = rawexp[(rawexp[var] < ul) & (rawexp[var] > ll)]
            # print('after', rawexp)
            rawptime = [float(i) for i in rawexp[var].tolist()]
            halfptime = min(rawptime) + ((max(rawptime) - min(rawptime)) / 2)
            start_index = []
            end_index = []
            for i in range(len(rawptime)):
                if rawptime[i] > halfptime:
                    end_index.append(i)
                else:
                    start_index.append(i)
            diff_results = []
            ps = []
            for gene in rawexp.columns.tolist()[9:]:
                if gene in ['psuedotime_mean', 'psuedotime_med', 'cell_ids']:
                    continue
                start_gene = []
                end_gene = []
                generaw = [float(i) for i in rawexp[gene].tolist()]
                for i in range(len(generaw)):
                    if i in end_index:
                        end_gene.append(generaw[i])
                    else:
                        start_gene.append(generaw[i])
                s = stats.ttest_ind(start_gene, end_gene)
                pval = s[1]
                teststat = s[0]
                diff = np.mean(end_gene) - np.mean(start_gene)
                if pval == pval:
                    diff_results.append([gene, diff, teststat, pval])
                    ps.append(pval)
            if len(ps) > 0:
                reject, fdr, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(ps, alpha=0.05,
                                                                                                 method='fdr_bh',
                                                                                                 is_sorted=False,
                                                                                                 returnsorted=False)
                for i in range(len(diff_results)):
                    diff_results[i].append(fdr[i])
                diff_results.sort(key=operator.itemgetter(4))
                folder = self.expression_path + '/' + '_'.join(branch.split(' -> ')) + '/'
                writeResults(diff_results, branch, 'raw_diffexp', folder)


########################################################

def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.

    Written by: Maulding
    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else:
        myCommandLine = CommandLine(
            myCommandLine)  # interpret the list passed from the caller of main as the commandline.

    try:

        start_time = time.time()
        pool = Pool(processes=myCommandLine.args.threads)

        all_dirs = []
        for d in os.listdir(myCommandLine.args.expression_path):
            if '.txt' not in d:
                if myCommandLine.args.expression_path[-1] == '/':
                    dir = myCommandLine.args.expression_path + d
                else:
                    dir = myCommandLine.args.expression_path + '/' + d
                all_dirs.append(dir)
        empty = [[] for _ in range(len(all_dirs))]

        result = pool.starmap(differential_genes, zip(all_dirs, empty))
        rawdiff(myCommandLine.args.expression_path, myCommandLine.args.json_dictionary)

        print("--- %s seconds ---" % (time.time() - start_time))
        if myCommandLine.args.expression_path[-1] == '/':
            f = open(myCommandLine.args.expression_path + 'diffexp_time.txt', 'x')
        else:
            f = open(myCommandLine.args.expression_path + '/' + 'diffexp_time.txt', 'x')
        f.write("--- %s seconds ---" % (time.time() - start_time))
        f.close()

    except Usage as err:
        print(err.msg)


if __name__ == "__main__":
    main()
