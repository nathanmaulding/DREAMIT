########### CommandLine ############################

class CommandLine():

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.

        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )
        self.parser.add_argument('-factor', help='input the factor you want to generate a plot from (e.g. MYC)')
        self.parser.add_argument('-method', default='All', help='input the method you want to generate a plot from (e.g. Pearson)')
        self.parser.add_argument('-branch', help='input the branch you are generating the plot from (e.g. 1_3)')
        self.parser.add_argument('-path', help='path to the branch eg /path/to/dreamit_output/branch_name/')

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

import numpy as np
import matplotlib
from numpy import nan
matplotlib.use('PS')
import matplotlib.pyplot as plt


#####################################################
## TO DO ##
# pick appropriate rolling examples

#####################################################

class figure:

    def __init__(self, branch, factor, method, path):
        self.branch = branch
        self.factor = factor
        self.method = method
        self.path = path

    def openFile(self, suffix):
        path = self.path + suffix
        data = open(path, 'r')
        figDict = ''
        for i, line in enumerate(data):
            line.split()
            if i >= 3:
                figDict = figDict + line
        figDict = eval(figDict)
        return figDict

    def plotting(self):
        ###################
        '''Needs to be adjusted to make non, pruned, and random plots'''
        ###################
        if self.method == 'Pearson' or self.method == 'Spearman' or self.method == 'All':
            nonrandom = self.openFile('Correlation_dictionary.txt')
            if self.factor not in nonrandom:
                print(self.factor, "is not in the dictionary. Double check spelling and letter casing.")
                return
            ####### unravel the data #########
            heat_arrs = nonrandom[self.factor][0]
            heat_y = nonrandom[self.factor][1]
            size = nonrandom[self.factor][2]

            ###### unravel pruned data #######
            pruned = self.openFile('Pruned_Pearson_Correlation_dictionary.txt')
            psize = size
            pheat_y = []
            pheat_arrs = []
            for i in range(len(heat_y)):
                if heat_y[i] in pruned[self.factor]:
                    pheat_y.append(heat_y[i])
                    pheat_arrs.append(heat_arrs[i])

            ##### unravel random data ########
            #random = self.openFile('Random_Correlation_dictionary.txt')
            rheat_arrs = nonrandom['Random1'][0] + nonrandom['Random2'][0] + nonrandom['Random3'][0]
            rheat_y = nonrandom['Random1'][1] + nonrandom['Random2'][1] + nonrandom['Random3'][1]
            rsize = nonrandom['Random1'][2]

            ##### Pearson heatmap #####
            def heatmapPNG(arr, y, size, branch, factor, type):
                fig, ax = plt.subplots()
                im = ax.imshow(arr)
                ax.set_yticks(np.arange(len(y)))
                ax.set_xticks(np.arange(size))
                ax.set_yticklabels(y)
                ax.set_xticklabels(range(1, size + 1))  # should be changed to label ticks with pseudotime
                ax.set_xlabel('bins')
                plt.colorbar(im)
                ax.set_title(branch + ' ' + factor + ' to targets')
                plt.savefig(branch + '_' + factor + '_' + type + '_heatmap.png')
                plt.clf()
                fig.clear()

            def expressionPlot(dictionary, branch, factor, type):
                if factor == 'Random':
                    for f in ['Random1', 'Random2', 'Random3']:
                        for i, gene_name in enumerate(dictionary[f][1]):
                            xplot = [x for x in range(1, dictionary[f][2] + 1)]  # should be changed to psuedotime
                            yplot = dictionary[f][0][i]
                            plt.plot(xplot, yplot, label=gene_name)
                else:
                    for i, gene_name in enumerate(dictionary[factor][1]):
                        xplot = [x for x in range(1, dictionary[factor][2] + 1)] # should be changed to psuedotime
                        yplot = dictionary[factor][0][i]
                        if gene_name == factor:
                            plt.plot(xplot, yplot, label=gene_name, color='black', linewidth=3)
                        else:
                            plt.plot(xplot, yplot, label=gene_name, alpha=0.5)
                plt.xlabel("bins")
                plt.ylabel("Normalized Expression")
                plt.legend()
                plt.title(branch + ' ' + factor + ' and target expression across bins')
                plt.savefig(branch + '_' + factor + '_' + type + '_expression.png')
                plt.clf()

            def prunedExpression(pruned, non, branch, factor, type):
                for i, gene_name in enumerate(non[factor][1]):
                    if gene_name in pruned[factor]:
                        xplot = [x for x in range(1, non[factor][2] + 1)] # should be changed to psuedotime
                        yplot = non[factor][0][i]
                        if gene_name == factor:
                            plt.plot(xplot, yplot, label=gene_name, color='black', linewidth=3)
                        else:
                            plt.plot(xplot, yplot, label=gene_name, alpha=0.5)
                plt.xlabel("bins")
                plt.ylabel("Normalized Expression")
                plt.legend()
                plt.title(branch + ' ' + factor + ' and target expression across bins')
                plt.savefig(branch + '_' + factor + '_' + type + '_expression.png')
                plt.clf()

            ###### Make heatmaps #######
            heatmapPNG(heat_arrs, heat_y, size, self.branch, self.factor, "Non_Random")
            heatmapPNG(pheat_arrs, pheat_y, psize, self.branch, self.factor, "Pruned")
            heatmapPNG(rheat_arrs, rheat_y, rsize, self.branch, self.factor, "Random")

            ##### Make expression plot #####
            expressionPlot(nonrandom, self.branch, self.factor, "Non_Random")
            expressionPlot(nonrandom, self.branch, 'Random', "Random")
            prunedExpression(pruned, nonrandom, self.branch, self.factor, "Pruned")


        if self.method == 'DTW' or self.method == 'All':
            nonrandom = self.openFile('DTW_dictionary.txt')
            if self.factor not in nonrandom:
                print(self.factor, "is not in the dictionary. Double check spelling and letter casing.")
                return

            ### get random data ####
            #random = self.openFile('Random_DTW_dictionary.txt')

            ### get pruned data ###
            pruned = self.openFile('Pruned_Dynamic_Time_Warping_dictionary.txt')

            ### unravel the data ###
            def unravelAndPlot(factor, dictionary, branch, type):
                plt.xlabel('targets')
                plt.ylabel(factor)
                if type != 'Random':
                    for target_i in dictionary[factor]:
                        target = target_i[0][0]
                        xpath = target_i[1]
                        ypath = target_i[2]
                        plt.plot(xpath, ypath, label=target)
                        plt.legend()
                if factor == 'Random':
                    for f in ['Random1', 'Random2', 'Random3']:
                        for target_i in dictionary[f]:
                            target = target_i[0][0]
                            xpath = target_i[1]
                            ypath = target_i[2]
                            plt.plot(xpath, ypath, label=target)
                plt.title(branch + ' ' + factor + type + ' DTW Paths')
                plt.savefig(branch + '_' + factor + type + '_dtwPaths.png')
                plt.clf()

            def prunedPlotting(factor, pruned_dictionary, nonrand_dictionary, branch, type):
                plt.xlabel('targets')
                plt.ylabel(factor)
                for target_i in nonrand_dictionary[factor]:
                    target = target_i[0][0]
                    xpath = target_i[1]
                    ypath = target_i[2]
                    if target in pruned_dictionary[factor]:
                        plt.plot(xpath, ypath, label=target)
                plt.legend()
                plt.title(branch + ' ' + factor + type + ' DTW Paths')
                plt.savefig(branch + '_' + factor + type + '_dtwPaths.png')
                plt.clf()

            ##### make image #####
            unravelAndPlot(self.factor, nonrandom, self.branch, 'Non_Random')
            unravelAndPlot('Random', nonrandom, self.branch, 'Random')
            prunedPlotting(self.factor, pruned, nonrandom, self.branch, 'Pruned')


        if self.method == 'MI' or self.method == 'All':
            nonrandom = self.openFile('Mutual_Information_dictionary.txt')
            if self.factor not in nonrandom:
                print(self.factor, "is not in the dictionary. Double check spelling and letter casing.")
                return
            #random = self.openFile('Random_Mutual_Information_dictionary.txt')
            pruned = self.openFile('Pruned_Mutual_Information_dictionary.txt')

            ### unravel the data ###
            def unravelAndPlot(factor, dictionary, branch, type):
                plt.xlabel('targets')
                plt.ylabel('Mutual Information Score')
                targets = []
                targetsMI = []
                if factor == 'Random':
                    for f in ['Random1', 'Random2', 'Random3']:
                        for target_i in dictionary[f]:
                            targets.append(target_i[0][0])
                            targetsMI.append(target_i[1])
                else:
                    for target_i in dictionary[factor]:
                        targets.append(target_i[0][0])
                        targetsMI.append(target_i[1])
                plt.bar(targets, targetsMI, color='black')
                plt.title(branch + ' ' + factor + type + ' Mutual Information')
                plt.savefig(branch + '_' + factor + type + '_MI.png')
                plt.clf()

            def prunedPlotting(factor, pruned_dictionary, nonrand_dictionary, branch, type):
                plt.xlabel('targets')
                plt.ylabel('Mutual Information Score')
                targets = []
                targetsMI = []
                for target_i in nonrand_dictionary[factor]:
                    if target_i[0][0] in pruned_dictionary[factor]:
                        targets.append(target_i[0][0])
                        targetsMI.append(target_i[1])
                plt.bar(targets, targetsMI, color='black')
                plt.title(branch + ' ' + factor + type + ' Mutual Information')
                plt.savefig(branch + '_' + factor + type + '_MI.png')
                plt.clf()

            ##### make image #####
            unravelAndPlot(self.factor, nonrandom, self.branch, 'Non_Random')
            unravelAndPlot('Random', nonrandom, self.branch, 'Random')
            prunedPlotting(self.factor, pruned, nonrandom, self.branch, 'Pruned')


        if self.method == 'Rolling' or self.method == 'All':

            def rollPlot(dictionary, factor, branch, type, last, ylab):
                fig, axs = plt.subplots(2)
                if factor == 'Random':
                    for f in ['Random1', 'Random2', 'Random3']:
                        for tar_i in dictionary[f]:
                            target = tar_i[0][0]
                            num_of_bins = tar_i[4]
                            mintime = tar_i[5]
                            maxtime = tar_i[6]
                            increment = (maxtime - mintime) / num_of_bins  # 35.4 / 20 = 1.77
                            x = tar_i[1]  # -36 to -4
                            y1 = tar_i[2]  #
                            y2 = tar_i[3]
                            xplot = []
                            xp = -(maxtime - mintime)
                            for i in range(len(x)):
                                # if (x[i] + num_of_bins) <= 0:
                                xp += increment
                                xplot.append(xp)
                            axs[0].plot(xplot, y1, label=target)
                            axs[1].plot(xplot, y2, label=target)
                else:
                    for tar_i in dictionary[factor]:
                        target = tar_i[0][0]
                        num_of_bins = tar_i[4]  # only plot the lagged pseudotime and subtract the number of bins from the x value
                        mintime = tar_i[5]
                        maxtime = tar_i[6]
                        increment = (maxtime - mintime) / num_of_bins # 35.4 / 20 = 1.77
                        x = tar_i[1]  # -36 to -4
                        y1 = tar_i[2]  #
                        y2 = tar_i[3]
                        xplot = []
                        xp = -(maxtime - mintime)
                        for i in range(len(x)):
                            # if (x[i] + num_of_bins) <= 0:
                            xp += increment
                            xplot.append(xp)
                        axs[0].plot(xplot, y1, label=target, alpha=0.5)  # ValueError: x and y must have same first dimension, but have shapes (14,) and (28,)
                        axs[1].plot(xplot, y2, label=target, alpha=0.5)
                ########### R^2 ( 2 plot lagged R^2) distribution plot #################
                axs[0].set(ylabel=ylab)
                axs[1].set(xlabel='Pseudotime Delay Before Relationship', ylabel='-log10(pvalue)')
                #axs[0].set_xlim(right=0)
                #axs[1].set_xlim(right=0)
                if type != 'Random':
                    plt.legend()
                axs[0].set(title=branch + ' ' + factor + ' ' + ' '.join((type + '_' + last).split('_')))
                plt.savefig(branch + '_' + factor + type + '_' + last + '.png')
                plt.clf()

            def prunedRolling(nonrandom_dict, pruned_dict, factor, branch, type, last, ylab):
                targets = []
                x = []
                y1 = []
                y2 = []
                for tar_i in nonrandom_dict[factor]:
                    if tar_i[0][0] in pruned_dict[factor]:
                        num_of_bins = int(tar_i[4])
                        mintime = tar_i[5]
                        maxtime = tar_i[6]
                        increment = (maxtime - mintime) / num_of_bins  # 35.4 / 20 = 1.77
                        targets.append(tar_i[0][0])
                        xp = []
                        xs = -(maxtime - mintime)
                        for i in tar_i[1]:
                            xs += increment
                            xp.append(xs)
                        x.append(xp)
                        y1.append(tar_i[2])
                        y2.append(tar_i[3])
                fig, axs = plt.subplots(2)
                for i in range(len(targets)):
                    axs[0].plot(x[i], y1[i], label=targets[i], alpha=0.5)
                    axs[1].plot(x[i], y2[i], label=targets[i], alpha=0.5)
                axs[0].set(ylabel=ylab)
                axs[1].set(xlabel='Pseudotime Delay Before Relationship', ylabel='-log10(pvalue)')
                #axs[0].set_xlim(right=0)
                #axs[1].set_xlim(right=0)
                plt.legend()
                axs[0].set(title=branch + ' ' + factor + ' ' + ' '.join((type + '_' + last).split('_')))
                plt.savefig(branch + '_' + factor + type + '_' + last + '.png')
                plt.clf()

            ########## Pearson #############
            ################################
            nonrandom = self.openFile('Rolling_Pearson_dictionary.txt')
            if self.factor not in nonrandom:
                print(self.factor, "is not in the dictionary. Double check spelling and letter casing.")
                return

            ##### get random data #####
            #random = self.openFile('Random_Rolling_Pearson_dictionary.txt')
            pruned = self.openFile('Pruned_Rolling_Pearson_Correlation_dictionary.txt')
            ##### make images #####
            rollPlot(nonrandom, self.factor, self.branch, 'Non_Random', 'Rolling_Pearson_Correlation', 'R^2')
            rollPlot(nonrandom, 'Random', self.branch, 'Random', 'Rolling_Pearson_Correlation', 'R^2')
            prunedRolling(nonrandom, pruned, self.factor, self.branch, 'Pruned', 'Rolling_Pearson_Correlation', 'R^2')


            ########## Spearman #############
            #################################
            nonrandom = self.openFile('Rolling_Spearman_dictionary.txt')
            if self.factor not in nonrandom:
                print(self.factor, "is not in the dictionary. Double check spelling and letter casing.")
                return

            ##### get random data #####
            #random = self.openFile('Random_Rolling_Spearman_dictionary.txt')
            pruned = self.openFile('Pruned_Rolling_Spearman_Correlation_dictionary.txt')
            ##### make images #####
            rollPlot(nonrandom, self.factor, self.branch, 'Non_Random', 'Rolling_Spearman_Correlation', 'S^2')
            rollPlot(nonrandom, 'Random', self.branch, 'Random', 'Rolling_Spearman_Correlation', 'S^2')
            prunedRolling(nonrandom, pruned, self.factor, self.branch, 'Pruned', 'Rolling_Spearman_Correlation', 'S^2')


            def rollPlotOther(dictionary, factor, branch, type, last, ylab):
                fig, axs = plt.subplots()
                if factor == 'Random':
                    for f in ['Random1', 'Random2', 'Random3']:
                        for tar_i in dictionary[f]:
                            target = tar_i[0][0]
                            num_of_bins = tar_i[3]
                            mintime = tar_i[4]
                            maxtime = tar_i[5]
                            increment = (maxtime - mintime) / num_of_bins
                            x = tar_i[1]  # -36 to -4
                            y1 = tar_i[2]
                            xplot = []
                            for i in range(len(x)):
                                # if (x[i] + num_of_bins) <= 0:
                                xpoint = x[i] + num_of_bins  # (-36 + 20) * 1.75 = -32
                                xplot.append(xpoint)
                            axs.plot(xplot, y1, label=target, alpha=0.5)  # TypeError: 'AxesSubplot' object is not subscriptable
                else:
                    for tar_i in dictionary[factor]:
                        target = tar_i[0][0]
                        num_of_bins = tar_i[3]  # only plot the lagged pseudotime and subtract the number of bins from the x value
                        mintime = tar_i[4]
                        maxtime = tar_i[5]
                        increment = (maxtime - mintime) / num_of_bins
                        x = tar_i[1]  # -36 to -4
                        y1 = tar_i[2]
                        xplot = []
                        for i in range(len(x)):
                            # if (x[i] + num_of_bins) <= 0:
                            xpoint = x[i] + num_of_bins  # (-36 + 20) * 1.75 = -32
                            xplot.append(xpoint)
                        axs.plot(xplot, y1, label=target, alpha=0.5)  # TypeError: 'AxesSubplot' object is not subscriptable
                ########### R^2 ( 2 plot lagged R^2) distribution plot #################
                axs.set(xlabel='bins moved', ylabel=ylab)
                if type != 'Random':
                    plt.legend()
                axs.set(title=branch + ' ' + factor + ' ' + ' '.join((type + '_' + last).split('_')))
                plt.savefig(branch + '_' + factor + type + '_' + last + '.png')
                plt.clf()

            def prunedRollingOther(nonrandom_dict, pruned_dict, factor, branch, type, last, ylab):
                targets = []
                x = []
                y1 = []
                y2 = []
                for tar_i in nonrandom_dict[factor]:
                    if tar_i[0][0] in pruned_dict[factor]:
                        num_of_bins = int(tar_i[3])
                        targets.append(tar_i[0][0])
                        xp = []
                        for i in tar_i[1]:
                            xp.append(int(i) + num_of_bins)
                        x.append(xp)
                        y1.append(tar_i[2])
                fig, axs = plt.subplots()
                for i in range(len(targets)):
                    axs.plot(x[i], y1[i], label=targets[i], alpha=0.5)
                axs.set(xlabel='bins moved', ylabel=ylab)
                plt.legend()
                axs.set(title=branch + ' ' + factor + ' ' + ' '.join((type + '_' + last).split('_')))
                plt.savefig(branch + '_' + factor + type + '_' + last + '.png')
                plt.clf()

            ########## Mutual Information #############
            ###########################################
            nonrandom = self.openFile('Rolling_MI_dictionary.txt')
            if self.factor not in nonrandom:
                print(self.factor, "is not in the dictionary. Double check spelling and letter casing.")
                return

            ##### get random data #####
            #random = self.openFile('Random_Rolling_MI_dictionary.txt')
            pruned = self.openFile('Pruned_Rolling_Mutual_Information_dictionary.txt')
            ##### make images #####
            rollPlotOther(nonrandom, self.factor, self.branch, 'Non_Random', 'Rolling_Mutual_Information', 'MI Content')
            rollPlotOther(nonrandom, 'Random', self.branch, 'Random', 'Rolling_Mutual_Information', 'MI Content')
            prunedRollingOther(nonrandom, pruned, self.factor, self.branch, 'Pruned', 'Rolling_Mutual_Information', 'MI Content')


            ########## DTW #############
            ############################
            nonrandom = self.openFile('Rolling_DTW_dictionary.txt')
            if self.factor not in nonrandom:
                print(self.factor, "is not in the dictionary. Double check spelling and letter casing.")
                return

            ##### get random data #####
            #random = self.openFile('Random_Rolling_DTW_dictionary.txt')
            pruned = self.openFile('Pruned_Rolling_DTW_dictionary.txt')
            ##### make images #####
            rollPlotOther(nonrandom, self.factor, self.branch, 'Non_Random', 'Rolling_DTW', 'DTW Cost')
            rollPlotOther(nonrandom, 'Random', self.branch, 'Random', 'Rolling_DTW', 'DTW Cost')
            prunedRollingOther(nonrandom, pruned, self.factor, self.branch, 'Pruned', 'Rolling_DTW', 'DTW Cost')


########################################################

def main(myCommandLine=None):

    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else:
        myCommandLine = CommandLine(
            myCommandLine)  # interpret the list passed from the caller of main as the commandline.

    try:

        result = figure(myCommandLine.args.branch, myCommandLine.args.factor, myCommandLine.args.method, myCommandLine.args.path).plotting()

    except Usage as err:
        print(err.msg)


if __name__ == "__main__":
    main()