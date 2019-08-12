#!/usr/bin/env python3
### CollectConsResults.py ---
##
## Filename: CollectConsResults.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Sep  6 09:02:18 2018 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Thu Jan 17 16:19:16 2019 (+0100)
##           By: Joerg Fallmann
##     Update #: 118
## URL:
## Doc URL:
## Keywords:
## Compatibility:
##
######################################################################
##
### Commentary:
###import os, sys, inspect
# # realpath() will make your script run, even if you symlink it :)
# cmd_folder = os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) ))
# if cmd_folder not in sys.path:
#     sys.path.insert(0, cmd_folder)
#
# # Use this if you want to include modules from a subfolder
# cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) )),"lib")
# if cmd_subfolder not in sys.path:
#     sys.path.insert(0, cmd_subfolder)
#
# # Info:
# # cmd_folder = os.path.dirname(os.path.abspath(__file__)) # DO NOT USE __file__ !!!
# # __file__ fails if the script is called in different ways on Windows.
# # __file__ fails if someone does os.chdir() before.
# # sys.argv[0] also fails, because it doesn't not always contains the path.
##
##
######################################################################
##
### Change Log:
##
##
######################################################################
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>.
##
######################################################################
##
### Code:

import os, sys, inspect
import glob
import argparse
from io import StringIO
import subprocess
import gzip
import importlib
import pprint
import traceback as tb
#numpy and matplolib and pyplot
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle # need that for histogram legend, don't ask
import numpy as np
#collections
from collections import Counter
from collections import defaultdict
import heapq
from operator import itemgetter
from natsort import natsorted, ns
#multiprocessing
import multiprocessing
from multiprocessing import Manager, Pool
#Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq
#load own modules
cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) )),"../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
from Collection import *

#sys.path=[str(os.getenv('HOME')+"/Work/Scripts/Python/")] + sys.path
#sys.path=[str(os.getenv('HOME')+"/Work/Scripts/Python/lib")] + sys.path
#sys.path=[str(os.getenv('HOME')+"/Work/Scripts/Python/Folding")] + sys.path

#parse args
def parseargs():
    parser = argparse.ArgumentParser(description='Calculate the regions with highest accessibility diff for given Sequence Pattern')
    parser.add_argument("-p", "--pattern", type=str, default='250,150', help='Pattern for files and window, e.g. Seq1_30,250')
    parser.add_argument("-c", "--cutoff", type=float, default=.2, help='Cutoff for the definition of pairedness, if set to e.g. 0.2 it will mark all regions with probability of being unpaired >= cutoff as unpaired')
    parser.add_argument("-b", "--border", type=str, default='', help='Cutoff for the minimum change between unconstraint and constraint structure, regions below this cutoff will not be returned as list of regions with most impact on structure. If not defined, will be calculated from folding the sequence of interest at temperature range 30-44.')
    parser.add_argument("-u", "--ulimit", type=int, default=1, help='Stretch of nucleotides used during plfold run (-u option)')
    parser.add_argument("-r", "--roi", type=str, default=None, help='Define Region of Interest that will be compared')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory to write to')
    parser.add_argument("-g", "--genes", type=str, help='Genomic coordinates bed for genes, either standard bed format or AnnotateBed.pl format')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processed to run this job with, only important of no border is given and we need to fold')
    return parser.parse_args()

def screen_genes(pat, cutoff, border, ulim, procs, roi, outdir, genes):
    try:
        #set path for output
        if outdir:
            printlog('Printing to ' + outdir)
            if not os.path.isabs(outdir):
                outdir =  os.path.abspath(outdir)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
        else:
            outdir = os.path.abspath(os.getcwd())

        pattern = pat.split(sep=',')
        window = int(pattern[0])
        span = int(pattern[1])

        genecoords = parse_annotation_bed(genes) #get genomic coords to print to bed later, should always be just one set of coords per gene

        # Create process pool with processes
        num_processes = procs or 1
        pool = multiprocessing.Pool(processes=num_processes, maxtasksperchild=1)

        for goi in genecoords:

            printlog('Working on ' + goi)
            gs, ge = map(int, genecoords[goi][0].split(sep='-'))

            #get files with specified pattern
            raw = os.path.abspath(os.path.join(goi, goi + '*_raw_*' + str(window) + '_' + str(span) + '.gz'))
            paired = os.path.abspath(os.path.join(goi, 'StruCons_' + goi + '*_diffnu_*' + str(window) + '_' + str(span) + '.gz'))
            unpaired = os.path.abspath(os.path.join(goi, 'StruCons_' + goi + '*_diffnp_*' + str(window) + '_' + str(span) + '.gz'))

            #search for files
            r = natsorted(glob.glob(raw), key=lambda y: y.lower())
            p = natsorted(glob.glob(paired), key=lambda y: y.lower())
            u = natsorted(glob.glob(unpaired), key=lambda y: y.lower())
            #c = natsorted(glob.glob(cutf), key=lambda y: y.lower())

            #get absolute path for files
            nocons = []

            raw = [os.path.abspath(i) for i in r]
            paired = [os.path.abspath(i) for i in p]
            unpaired = [os.path.abspath(i) for i in u]

            for i in range(len(r)):
                pool.apply_async(judge_diff, args=(raw[i], u[i], p[i], gs, ge, ulim, cutoff, border, outdir),error_callback=eprint)

        pool.close()
        pool.join()

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)


def judge_diff(raw, u, p, gs, ge, ulim, cutoff, border, outdir):

    try:
        goi, chrom, strand, cons, reg, f, window, span = map(str,os.path.basename(raw).split(sep='_'))
        cs, ce = map(int, cons.split(sep='-'))
        ws, we = map(int, reg.split(sep='-'))

        cs = cs - ws #fit to window and make 0-based
        ce = ce - ws #fit to window and make 0-based
        if strand is not '-':
            ws = ws + gs -1 #get genomic coords
            we = we + gs
        else:
            wst = ws            # temp ws for we calc
            ws = ge - we - 1#get genomic coords
            we = ge - wst

        border1, border2 = map(float,border.split(',')) #defines how big a diff has to be to be of importance

        printlog('Continuing calculation with cutoff: ' + str(cutoff) + ' and borders ' + str(border1) + ' and ' + str(border2))

        out = {}
        out['p'] = []
        out['u'] = []

        noc = pl_to_array(raw, ulim-1)

        if abs(noc[ce]) > cutoff:
            uc = pl_to_array(u, ulim-1)
            pc = pl_to_array(p, ulim-1)

            for pos in range(len(uc)):
                if pos not in range(cs,ce+1):
                    if border1 < uc[pos] and uc[pos] < border2:
                        if ce < pos:# get distance up or downstream
                            dist = (pos - ce) * -1 # no -1 or we have 0 overlap
                        else:
                            dist = cs - pos

                        if strand is not '-':
                            gpos = pos + ws
                            gend = gpos + ulim
                            gcons = str(cs+ws)+'-'+str(ce+ws)
                        else:
                            gpos = we - pos - ulim-1
                            gend = gpos + ulim
                            gcons = str(we-ce-1)+'-'+str(we-cs)

#                        gcons = str(cs+ws)+'-'+str(ce+ws)

                        out['u'].append('\t'.join([str(chrom), str(gpos), str(gend), str(goi) + '|' + str(cons) + '|' + str(gcons), str(uc[pos]), str(strand), str(dist), str(noc[pos])]))

                    if border1 < pc[pos] and pc[pos] < border2:
                        if ce < pos:# get distance up or downstream
                            dist = (pos - ce) * -1 # no -1 or we have 0 overlap
                        else:
                            dist = cs - pos

                        if strand is not '-':
                            gpos = pos + ws
                            gend = gpos + ulim
                            gcons = str(cs+ws)+'-'+str(ce+ws)
                        else:
                            gpos = we - pos - ulim-1
                            gend = gpos + ulim
                            gcons = str(we-ce-1)+'-'+str(we-cs)

                        out['p'].append('\t'.join([str(chrom), str(gpos), str(gend), str(goi) + '|' + str(cons) + '|' + str(gcons), str(pc[pos]), str(strand), str(dist), str(noc[pos])]))

        savelists(out, outdir)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def savelists(out, outdir):
    try:
        if len(out['u']) > 0:
            o = gzip.open(os.path.abspath(os.path.join(outdir, 'Collection_unpaired.bed.gz')), 'ab')
            o.write(bytes('\n'.join(out['u']),encoding='UTF-8'))
            o.write(bytes('\n',encoding='UTF-8'))
        if len(out['p']) > 0:
            o = gzip.open(os.path.abspath(os.path.join(outdir, 'Collection_paired.bed.gz')), 'ab')
            o.write(bytes('\n'.join(out['p']),encoding='UTF-8'))
            o.write(bytes('\n',encoding='UTF-8'))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

####################
####    MAIN    ####
####################
if __name__ == '__main__':
    args=parseargs()
    screen_genes(args.pattern, args.cutoff, args.border, args.ulimit, args.procs, args.roi, args.outdir, args.genes)
######################################################################
### CollectConsResults.py ends here
