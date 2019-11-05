#!/usr/bin/env python3
## CollectWindowResults.py ---
##
## Filename: CollectWindowResults.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Aug 15 13:49:46 2019 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Tue Nov  5 13:59:53 2019 (+0100)
##           By: Joerg Fallmann
##     Update #: 62
## URL:
## Doc URL:
## Keywords:
## Compatibility:
##
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
### IMPORTS
import os, sys, inspect

##load own modules
from lib.Collection import *
from lib.logger import makelogdir, setup_multiprocess_logger
# Create log dir
makelogdir('logs')
# Define loggers
scriptname=os.path.basename(__file__)
global streamlog, log           # global to ensure that later manipulation of loglevel is possible
streamlog = setup_multiprocess_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='WARNING')
log = setup_multiprocess_logger(name=scriptname, log_file='logs/'+scriptname, logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='WARNING')

##other modules
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

#parse args
def parseargs():
    parser = argparse.ArgumentParser(description='Calculate the regions with highest accessibility diff for given Sequence Pattern')
    parser.add_argument("-p", "--pattern", type=str, default='250,150', help='Pattern for files and window, e.g. Seq1_30,250')
    parser.add_argument("-b", "--border", type=str, default='', help='Cutoff for the minimum change between unconstraint and constraint structure, regions below this cutoff will not be returned as list of regions with most impact on structure. If not defined, will be calculated from folding the sequence of interest at temperature range 30-44.')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory to write to')
    parser.add_argument("-g", "--genes", type=str, help='Genomic coordinates bed for genes, either standard bed format or AnnotateBed.pl format')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processed to run this job with, only important of no border is given and we need to fold')
    parser.add_argument("--loglevel", type=str, default='WARNING', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def screen_genes(pat, border, procs, outdir, genes):

    logid = scriptname+'.screen_genes: '
    try:
        #set path for output
        if outdir:
            log.info(logid+'Printing to ' + outdir)
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

            log.info(logid+'Working on ' + goi)
            gs, ge = map(int, genecoords[goi][0].split(sep='-'))

            #get files with specified pattern
            paired = os.path.abspath(os.path.join(goi, goi + '*_pairedconstraint_*' + str(window) + '_' + str(span) + '.gz'))

            #search for files
            p = natsorted(glob.glob(paired), key=lambda y: y.lower())
            log.debug(logid+'Files found: ' + str(p))

            #get absolute path for files
            nocons = []

            paired = [os.path.abspath(i) for i in p]

            if not paired:
                log.warning(logid+'No output for gene '+str(goi)+' found, will skip!')
                continue

            try:
                for i in range(len(p)):
                    log.debug(logid+'Calculating file ' + str(p[i]))
                    pool.apply_async(calc, args=(p[i], gs, ge, border, outdir))
            except Exception as err:
                exc_type, exc_value, exc_tb = sys.exc_info()
                tbe = tb.TracebackException(
                    exc_type, exc_value, exc_tb,
                )
                log.error(logid+''.join(tbe.format()))

        pool.close()
        pool.join()

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def calc(p, gs, ge, border, outdir):

    logid = scriptname+'.calc_ddg: '
    try:
        goi, chrom, strand, cons, reg, window, span = map(str,os.path.basename(p).split(sep='_'))
        border1, border2 = map(float,border.split(',')) #defines how big a diff has to be to be of importance

        log.info(logid+'Continuing calculation with borders: ' + str(border1) + ' and ' + str(border2))

        out = defaultdict()
        ddgs = get_ddg(p)
        log.debug(logid+str(ddgs))

        RT = (-1.9872041*10**(-3))*(37+273.15)
        log.debug(logid+'RT is '+str(RT))

        for cons in ddgs:
            if not cons in out:
                out[cons] = list()
            ddg=calc_ddg(ddgs[cons])
            if ddg is not None:
                if ddg > border1 and ddg < border2:
                    dkd = math.exp(ddg/RT)
                    out[cons].append('\t'.join([str(chrom), str(gs), str(ge),  str(goi), str(ddg), str(strand), str(cons),str(dkd)+'\n']))
        if out:
            write_out(out, outdir)
        else:
            log.warning(logid+'No ddg above cutoffs for gene '+str(goi))
        return

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def write_out(out, outdir):

    logid = scriptname+'.savelist: '
    try:
        for cons in out:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            if not os.path.exists(os.path.abspath(os.path.join(outdir, 'Collection_window.bed.gz'))):
                with gzip.open(os.path.abspath(os.path.join(outdir, 'Collection_window.bed.gz')), 'wb') as o:
                    o.write(bytes('\n'.join(out[cons]),encoding='UTF-8'))
            else:
                with gzip.open(os.path.abspath(os.path.join(outdir, 'Collection_window.bed.gz')), 'ab') as o:
                    o.write(bytes('\n'.join(out[cons]),encoding='UTF-8'))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        if args.loglevel != 'WARNING':
          streamlog = setup_multiprocess_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
          log = setup_multiprocess_logger(name=scriptname, log_file='logs/'+scriptname, logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)

        log.info(logid+'Running '+scriptname+' on '+str(args.procs)+' cores')
        screen_genes(args.pattern, args.border, args.procs, args.outdir, args.genes)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

######################################################################
#
# CollectWindowResults.py ends here
