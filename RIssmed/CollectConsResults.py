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
## Last-Updated: Thu Feb  6 15:06:03 2020 (+0100)
##           By: Joerg Fallmann
##     Update #: 371
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
import math
from scipy.stats import zscore as zsc

def parseargs():
    parser = argparse.ArgumentParser(description='Calculate the regions with highest accessibility diff for given Sequence Pattern')
    parser.add_argument("-p", "--pattern", type=str, default='250,150', help='Pattern for files and window, e.g. Seq1_30,250')
    parser.add_argument("-c", "--cutoff", type=float, default=.2, help='Cutoff for the definition of pairedness, if set to e.g. 0.2 it will mark all regions with probability of being unpaired >= cutoff as unpaired')
    parser.add_argument("-b", "--border", type=str, default='', help='Cutoff for the minimum change between unconstraint and constraint structure, regions below this cutoff will not be returned as list of regions with most impact on structure. If not defined, will be calculated from folding the sequence of interest at temperature range 30-44.')
    parser.add_argument("-u", "--ulimit", type=int, default=1, help='Stretch of nucleotides used during plfold run (-u option)')
    parser.add_argument("-r", "--roi", type=str, default=None, help='Define Region of Interest that will be compared')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory to write to')
    parser.add_argument("-d", "--dir", type=str, default='', help='Directory to read from')
    parser.add_argument("-g", "--genes", type=str, help='Genomic coordinates bed for genes, either standard bed format or AnnotateBed.pl format')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processes to run this job with, only important of no border is given and we need to fold')
    parser.add_argument("--loglevel", type=str, default='WARNING', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")
    parser.add_argument("-w", "--padding", type=int, default=1, help='Padding around constraint that will be excluded from report, default is 1, so directly overlapping effects will be ignored')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def screen_genes(pat, cutoff, border, ulim, procs, roi, outdir, genes, padding):

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

        log.debug(logid+str(genecoords))

        # Create process pool with processes
        num_processes = procs or 1
        pool = multiprocessing.Pool(processes=num_processes, maxtasksperchild=1)

        for goi in genecoords:
            log.info(logid+'Working on ' + goi)
            gs, ge, gstrand = get_location(genecoords[goi][0])

            #get files with specified pattern
            raw = os.path.abspath(os.path.join(goi, goi + '*_raw_*' + str(window) + '_' + str(span) + '.npy'))
            unpaired = os.path.abspath(os.path.join(goi, 'StruCons_' + goi + '*_diffnu_*' + str(window) + '_' + str(span) + '.npy'))
            paired = os.path.abspath(os.path.join(goi, 'StruCons_' + goi + '*_diffnp_*' + str(window) + '_' + str(span) + '.npy'))

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

            if not raw or not paired or not unpaired:
                log.warning(logid+'Could not find files for Gene '+str(goi)+' and window '+str(window)+' and span '+str(span)+' Will skip')
                continue

            try:
                for i in range(len(r)):
                    pool.apply_async(judge_diff, args=(raw[i], u[i], p[i], gs, ge, gstrand, ulim, cutoff, border, outdir, padding))
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

def judge_diff(raw, u, p, gs, ge, gstrand, ulim, cutoff, border, outdir, padding):

    logid = scriptname+'.judge_diff: '
    try:
        goi, chrom, strand, cons, reg, f, window, span = map(str,os.path.basename(raw).split(sep='_'))
        span = span.split(sep='.')[0]
        cs, ce = map(int, cons.split(sep='-'))
        ws, we = map(int, reg.split(sep='-'))

        cs = cs - ws #fit to window and make 0-based
        ce = ce - ws #fit to window and make 0-based closed

        if 0 > any([cs,ce,ws,we]):
            raise Exception('One of '+str([cs,ce,ws,we])+ ' lower than 0! this should not happen for '+','.join([goi, chrom, strand, cons, reg, f, window, span]))

        if gstrand is not '-':
            ws = ws + gs - 2 #get genomic coords 0 based closed, ws and gs are 1 based
            we = we + gs - 2

        else:
            wst = ws     #temp ws for we calc
            ws = ge - we #get genomic coords 0 based closed, ge and we are 1 based
            we = ge - wst

        log.debug(logid+'DiffCoords: '+' '.join(map(str,[goi, chrom, strand, cons, reg, f, window, span, gs, ge, cs, ce, ws, we])))

        border1, border2 = map(float,border.split(',')) #defines how big a diff has to be to be of importance

        log.info(logid+'Continuing '+str(goi)+' calculation with cutoff: ' + str(cutoff) + ' and borders ' + str(border1) + ' and ' + str(border2))

        out = {}
        out['p'] = []
        out['u'] = []

        RT = (-1.9872041*10**(-3))*(37+273.15)
        log.debug(logid+'RT is '+str(RT))

        noc = pl_to_array(raw, ulim)

        mult = int((len(noc)/int(window))/2)
        log.debug(logid+'Multiplyer: '+str(mult))
        cws = int(window)*(mult-1)
        if cws < window:


        cwe = int(window)*(mult+1)+ulim-1
        if cwe > len(noc):
            cwe = len(noc)
        conswindow = (cws,cwe) #0-based half open
        log.debug(logid+'Constraint Window: '+str(conswindow))

        if abs(noc[ce]) > cutoff:
            uc = pl_to_array(u, ulim)  # This is the diffacc for unpaired constraint
            pc = pl_to_array(p, ulim)  # This is the diffacc for paired constraint

            log.debug(logid+'unpaired: '+str(u)+' and paired: '+str(p)+' Content: '+str(uc[ulim:ulim+10])+' test '+str(np.all(uc[ulim:ulim+10])))

            epsilon = 10 ** -50
            preaccu = noc+uc+epsilon
            preaccp = noc+pc+epsilon

            np.seterr(divide = 'ignore') #ignore 0 for logs
            nrgdiffu = np.array(RT*np.log(abs(uc)))
            nrgdiffp = np.array(RT*np.log(abs(pc)))
            np.seterr(divide = 'warn')

            #replace -inf with nan
            nrgdiffu[np.isneginf(nrgdiffu)] = np.nan
            nrgdiffp[np.isneginf(nrgdiffp)] = np.nan

            kdu = np.exp(nrgdiffu/RT)#math.exp(np.array(nrgdiffu//RT))) ###THIS IS BASICALLY ACCESSIBILITY AGAIN
            kdp = np.exp(nrgdiffp/RT)#math.exp(np.array(nrgdiffp//RT))) ###THIS IS BASICALLY ACCESSIBILITY AGAIN

            log.debug(logid+'NRG: '+str(nrgdiffu[:10]))
            log.debug(logid+'KD: '+str(kdu[:10])+' mean: '+str(np.nanmean(kdu))+' std: '+str(np.nanstd(kdu)))

            np.seterr(divide = 'ignore') #ignore 0 for logs
            zscoresu = np.array(np.divide(kdu - np.nanmean(kdu), np.nanstd(kdu,ddof=0), out=np.zeros_like(kdu - np.nanmean(kdu)), where=np.nanstd(kdu,ddof=0)!=0)) #np.array(zsc(kdu[~np.isnan(kdu)]))
            zscoresp = np.array(np.divide(kdp - np.nanmean(kdp), np.nanstd(kdp,ddof=0), out=np.zeros_like(kdp - np.nanmean(kdp)), where=np.nanstd(kdp,ddof=0)!=0)) #np.array((kdp - np.nanmean(kdp))/np.nanstd(kdp,ddof=0))#np.array(zsc(kdp[~np.isnan(kdp)]))
            np.seterr(divide = 'warn')

            #replace -inf with nan
            zscoresu[np.isneginf(zscoresu)] = np.nan
            zscoresp[np.isneginf(zscoresp)] = np.nan

            log.debug(logid+'zscore: '+str(zscoresu[:10]))

            '''
            Collect positions of interest with padding around constraint
            Constraints are influencing close by positions strongest so strong influence of binding there is expected
            '''

            log.debug(logid+'WINDOWS: '+str.join(' ',map(str,[goi,conswindow[0],conswindow[1]+1,strand,ws,cs,ce,we,str(cs+ws-1)+'-'+str(ce+ws),str(we-ce-1)+'-'+str(we-cs)])))

            for pos in range(conswindow[0],conswindow[1]+1):
                if pos not in range(cs-padding-ulim,ce+1+padding+ulim):
                    if border1 < uc[pos] and uc[pos] < border2:
                        if ce < pos:# get distance up or downstream
                            dist = (pos - ce) * -1 # no -1 or we have 0 overlap
                        else:
                            dist = cs - pos

                        if strand is not '-':
                            gpos = pos + ws #already 0-based
                            gend = gpos + ulim #0-based half-open
                            gcons = str(cs+ws)+'-'+str(ce+ws+1)
                        else:
                            gpos = we - pos - ulim #already 0-based
                            gend = gpos + ulim + 1 #0-based half-open
                            gcons = str(we-ce-1)+'-'+str(we-cs)

                        preacc = preaccu[pos] - epsilon
                        nrgdiff = nrgdiffu[pos]
                        kd = kdu[pos]
                        zscore = zscoresu[pos]

                        if not any([x is np.nan for x in [preacc,nrgdiff,kd,zscore]]):
                            out['u'].append('\t'.join([str(chrom), str(gpos), str(gend), str(goi) + '|' + str(cons) + '|' + str(gcons), str(uc[pos]), str(strand), str(dist), str(noc[pos]), str(preacc), str(nrgdiff), str(kd), str(zscore)]))

                    if border1 < pc[pos] and pc[pos] < border2:
                        if ce < pos:# get distance up or downstream
                            dist = (pos - ce) * -1 # no -1 or we have 0 overlap
                        else:
                            dist = cs - pos

                        if strand is not '-':
                            gpos = pos + ws #already 0-based
                            gend = gpos + ulim #0-based half-open
                            gcons = str(cs+ws)+'-'+str(ce+ws+1)
                        else:
                            gpos = we - pos - ulim #already 0-based
                            gend = gpos + ulim + 1 #0-based half-open
                            gcons = str(we-ce-1)+'-'+str(we-cs)

                        preacc = preaccp[pos] - epsilon
                        nrgdiff = nrgdiffp[pos]
                        kd = kdp[pos]
                        zscore = zscoresp[pos]

                        if not any([x is np.nan for x in [preacc,nrgdiff,kd,zscore]]):
                            out['p'].append('\t'.join([str(chrom), str(gpos), str(gend), str(goi) + '|' + str(cons) + '|' + str(gcons), str(pc[pos]), str(strand), str(dist), str(noc[pos]), str(preacc), str(nrgdiff), str(kd)]))

        savelists(out, outdir)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def savelists(out, outdir):

    logid = scriptname+'.savelist: '
    try:
        #if not os.path.isfile(os.path.abspath(os.path.join(outdir, 'Collection_unpaired.bed.gz'))):
        #    with gzip.open(os.path.abspath(os.path.join(outdir, 'Collection_unpaired.bed.gz')), 'ab') as o:
        #        o.write(bytes('\n'.join(out['u']),encoding='UTF-8'))
        #        o.write(bytes('\n',encoding='UTF-8'))
        #
        #if not os.path.isfile(os.path.abspath(os.path.join(outdir, 'Collection_paired.bed.gz'))):
        #    with gzip.open(os.path.abspath(os.path.join(outdir, 'Collection_paired.bed.gz')), 'ab') as o:
        #        o.write(bytes('\n'.join(out['p']),encoding='UTF-8'))
        #        o.write(bytes('\n',encoding='UTF-8'))
        #[str(chrom), str(gpos), str(gend), str(goi) + '|' + str(cons) + '|' + str(gcons), str(uc[pos]), str(strand), str(dist), str(noc[pos]), str(accdiff), str(nrgdiff), str(kd), str(zscore)]
        if len(out['u']) > 0:
            with gzip.open(os.path.abspath(os.path.join(outdir, 'Collection_unpaired.bed.gz')), 'ab') as o:
                o.write(bytes('\n'.join(out['u']),encoding='UTF-8'))
                o.write(bytes('\n',encoding='UTF-8'))
        if len(out['p']) > 0:
            with gzip.open(os.path.abspath(os.path.join(outdir, 'Collection_paired.bed.gz')), 'ab') as o:
                o.write(bytes('\n'.join(out['p']),encoding='UTF-8'))
                o.write(bytes('\n',encoding='UTF-8'))
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
            #          streamlog = setup_multiprocess_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
            log = setup_multiprocess_logger(name=scriptname, log_file='logs/'+scriptname, logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
        if args.dir != '':
            args.genes = os.path.abspath(args.genes)
            os.chdir(os.path.abspath(args.dir))

        log.info(logid+'Running '+scriptname+' on '+str(args.procs)+' cores')
        screen_genes(args.pattern, args.cutoff, args.border, args.ulimit, args.procs, args.roi, args.outdir, args.genes, args.padding)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

######################################################################
### CollectConsResults.py ends here
