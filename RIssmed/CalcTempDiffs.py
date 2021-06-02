#!/usr/bin/env python3
# CalcTempDiffs.py ---
#
# Filename: CalcConsDiffs.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Nov 21 12:43:42 2017 (+0100)
# Version:
# Package-Requires: ()
# Last-Updated: Mon Sep 24 13:01:21 2018 (+0200)
#           By: Joerg Fallmann
#     Update #: 721
# URL:
# Doc URL:
# Keywords:
# Compatibility:
#
#

# Commentary:
#
#
#
#

# Change Log:
#
#
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>.
#
#

# Code:
import glob
import gzip
from matplotlib import pyplot as plt
from collections import defaultdict
import multiprocessing
from multiprocessing import Manager
from natsort import natsorted
#Biopython stuff
from Bio import SeqIO

from RIssmed.RNAtweaks.RIssmedArgparsers import *
#from ConstraintPLfold.py import

#parse args
def parseargs():
    parser = argparse.ArgumentParser(description='Calculate the regions with highest accessibility diff for given Sequence Pattern')
    parser.add_argument("-p", "--pattern", type=str, default='*,*.svg', help='Pattern for files and window, e.g. Seq1_30,250')
    parser.add_argument("-c", "--cutoff", type=int, default=.4, help='Cutoff for the definition of pairedness, if set to e.g. 0.2 it will mark all regions with probability of being unpaired > cutoff as paired')
    parser.add_argument("-u", "--ulimit", type=int, default=1, help='Stretch of nucleotides used during plfold run (-u option)')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processed to run this job with')
    parser.add_argument("-s", "--sequence", type=str, default='', help='Folded Sequence in FASTA format')
    parser.add_argument("--plot", type=int, default=1, help='Plot the distribution of probability differences at each nucleotide over the given temperature range')
    return parser.parse_args()

def calctdiff(pat, cutoff, ulim, sequence, plot=None, procs=None):

    pattern = pat.split(sep=',')
    window = int(pattern[1])
#get files with specified pattern
    raw = pattern[0] + '*_raw_' + pattern[1] + '.gz'
    temp = 'TempCons_' + pattern[0] + '*_temp_'+pattern[1] + '.gz'
#search for files
    r = natsorted(glob.glob(raw), key=lambda y: y.lower())
    t = natsorted(glob.glob(temp), key=lambda y: y.lower())
#37 degree is default, we use this as raw in case that probtemp folding was used
    for f in t:
        if '37_temp' in f:
            p = f
            t.remove(f)
#get absolute path for files
#we start with the temp constraint files
    temp = [os.path.abspath(i) for i in t]
    tconstest = toarray(temp[0], ulim)

    if r: #this means that there is a raw file
        raw = os.path.abspath(r[0])
    elif p: #this means we take the 37 file as raw
        raw = os.path.abspath(p)
#Now we check if nocons and tcons have the same length, if not most likely from probcons folding, then we use the 37 degree folding
    else:
        print('There is no raw file and no 37 folding temp file to compensate, this will fail')
        return 0

    nocons = toarray(raw, ulim)

#get Sequences matching pattern
    nucs = []
    seq = ''

    if (os.path.isfile(sequence)):
        if '.gz' in sequence:
            seq = gzip.open(sequence,'rt')
        else:
            seq = open(sequence,'rt')
        for fa in SeqIO.parse(seq,'fasta'):
            if pattern[0] in fa.id:
                nucs = fa.seq
    else:
        print('No sequence information available, will not add 2nd x axis to plot')

#read in the files, define pairedness and calculate statistics for distance to constraint vs diff in pairedness

#now we go through all the constraints, get the constraint region and calculate the difference
#create a nested dict for over all max/min regions
#   maxdiffcons = defaultdict(
#       lambda: defaultdict(
#           lambda: defaultdict(list)
#           )
#       )

###### Multicore processing
    manager = Manager()
    difftcons = manager.list()
    trange = manager.list()
    num_processes = procs or 1
    pool = multiprocessing.Pool(processes=num_processes)
    for i in range(len(temp)):
        pool.apply_async(calcdiff, args=(temp[i], nocons, trange, difftcons, ulim))
    pool.close()
    pool.join()

#Prepare for plots and analysis
    difftcons = np.array(difftcons)
    mint = getlowest_list(trange,1)[0]
    maxt = gethighest_list(trange,1)[0]
    trange=[]
#Transpose rows to columns for boxplot and analysis later
    difftcons = difftcons.transpose()
    meanofmeans, meanofstds = analyze(difftcons)
    with open('Cutoffs_' + pattern[0] + '_' + pattern[1] + '_' + mint + '_' + maxt + '.txt', 'w') as handle:
        print('MEAN:\t' + str(meanofmeans) + '\n' + 'STD:\t' + str(meanofstds), file=handle)

    if plot:
        plottempfig(nocons, difftcons, nucs, pattern[0], pattern[1], mint, maxt)

    return meanofmeans, meanofstds

###### SUBS
def calcdiff(file, nocons, trange, difftcons, ulim):
    trange.append(file.split(sep='/')[-1].split(sep='_')[4])
    tcons = toarray(file, ulim)
# calc diff between raw and constraint
    td = nocons - tcons
    difftcons.append(td)

def analyze(difftcons):
    means = difftcons.mean(axis=0)
    stddevs = difftcons.std(axis=0)
    meanofmeans = means.mean()
    meanofstds  = stddevs.mean()
    return meanofmeans, meanofstds

def plottempfig(nocons, difftcons, nucs, p1, p2, mint, maxt):
### Plot this distribution
    width = 16/100*len(nocons)
    height = 9
    fig = plt.figure(figsize=(width, height),dpi=100)
    ax = plt.subplot(111)
    ax.plot(range(1,len(nocons)+1), nocons, 'b-')
#Add boxplot of changes
    for i in range(len(difftcons)):
        ax.boxplot(difftcons[i],positions = [i+1], notch = True, meanline=True, showmeans=True, showcaps=True, showbox=True, showfliers=False)
    ax.set_xlim(0, len(nocons)+1)
    ax.set_xticks(range(0,len(nocons)+1))
    ax.set_xticklabels(range(0,len(nocons)+1), rotation=45, ha="center", fontsize=6)
    ax2 = ax.twiny()
    if nucs:
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(range(0,len(nocons)+1))
        ax2.set_xticklabels((' '+nucs), ha="center", fontsize=6)
#Add ticks and labels
    plt.xlabel("Nucleotide", fontsize=16)
    plt.ylabel("Prob of being unpaired and changes", fontsize=16)
    fig.savefig('DiffTemp_' + p1 + '_' + p2 + '_' + mint + '_' + maxt + '.svg')

def getdist(a, conss, conse):
    dist = []
    for i in range(len(a)):
        if a[i] != 0:
            if (int(conss) - i) > 0:
                dist.append(min(int(conss) - i, int(conse) - i) + 1)
            else:
                dist.append(max(int(conss) - i, int(conse) - i) + 1)
    return dist

def probbin(a, cutoff):
    b = a
    for i in np.nditer(b, op_flags=['readwrite'], casting='no', op_dtypes=float):
        if i[...] >= cutoff:
            i[...] = 1
        else:
            i[...] = 0
    return b

def borderchange(raw, constraint, border):
    bigdiff = {}
    diff = raw - constraint
    for i in range(len(diff)):
        if diff[i] <= border*(-1) or diff[i] >= border:
            bigdiff[i] = diff[i]
    return bigdiff

def maxchange(raw, constraint):
    bigdiff = []
    diff = raw - constraint
    for i in range(len(diff)):
        if diff[i] <= border*(-1) or diff[i] >= border:
            txt = str(i)+'\t'+str(diff[i])
            bigdiff.append(txt)
    return bigdiff

def savelist(gothrough, what, const):
    out = []
    for x,y in gothrough.items():
        out.append(str(const+'\t'+str(x)+'\t'+str(y)))
    o = gzip.open('Impact_of_constraints_on_'+what+'.gz', 'ab')
    o.write(bytes('\n'.join(out),encoding='UTF-8'))

def savemax(gothrough, which, what):
    tempdict = defaultdict(list)
    templist = []
    for z in ['lowest', 'highest']:
        for x,y in gothrough[what][z].items():
            for a,b in y.items():
                templist.append(b)
                tempdict[b].append(str(str(x)+'\t'+str(a)))
        call = {getlowest_list:'lowest',gethighest_list:'highest'}
        printlist = []
        for key,val in call.items():
            if val == z:
                printlist.extend(key(templist, 100))

        out = []
        for x in printlist:
            out.append(str(''.join(tempdict[x])+'\t'+str(x)))
        o = gzip.open('Impact'+z+'_'+what+'.gz', 'wb+')
        o.write(bytes('\n'.join(out),encoding='UTF-8'))

####################
####    MAIN    ####
####################
if __name__ == '__main__':
    args=parseargs()
    calctdiff(args.pattern, args.cutoff, args.ulimit, args.sequence, args.plot, args.procs)

#
# CalcTempDiffs.py ends here
