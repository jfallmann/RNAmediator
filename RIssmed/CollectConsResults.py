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
## Last-Updated: Mon May 11 15:27:01 2020 (+0200)
##           By: Joerg Fallmann
##     Update #: 450
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
import os
import sys
import inspect
import glob
import gzip
import traceback as tb
import argparse
# multiprocessing
import multiprocessing
from multiprocessing import get_context
from multiprocessing import set_start_method
# numpy
import numpy as np
import shlex
# others
from natsort import natsorted
# Logging
import datetime
import logging
from lib.logger import makelogdir, makelogfile, listener_process, listener_configurer, worker_configurer
# load own modules
from lib.Collection import *
from lib.FileProcessor import *
from lib.RNAtweaks import *
from lib.NPtweaks import *

log = logging.getLogger(__name__)  # use module name
scriptname = os.path.basename(__file__).replace('.py', '')


def screen_genes(queue, configurer, level, pat, cutoff, border, ulim, procs, roi, outdir, dir, genes, padding):

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
            raw = os.path.abspath(os.path.join(dir, goi, goi + '*_raw_*' + str(window) + '_' + str(span) + '.npy'))
            unpaired = os.path.abspath(os.path.join(dir, goi, 'StruCons_' + goi + '*_diffnu_*' + str(window) + '_' + str(span) + '.npy'))
            paired = os.path.abspath(os.path.join(dir, goi, 'StruCons_' + goi + '*_diffnp_*' + str(window) + '_' + str(span) + '.npy'))

            log.debug(logid+'PATHS: '+str(raw)+'\t'+str(paired)+'\t'+str(unpaired))

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

            log.debug(logid+'PATHS: '+str(raw)+'\t'+str(paired)+'\t'+str(unpaired))

            if not raw or not paired or not unpaired:
                log.warning(logid+'Could not find files for Gene '+str(goi)+' and window '+str(window)+' and span '+str(span)+' Will skip')
                continue

            try:
                for i in range(len(r)):
                    pool.apply_async(judge_diff, args=(raw[i], u[i], p[i], gs, ge, gstrand, ulim, cutoff, border, outdir, padding), kwds={'queue':queue, 'configurer':configurer, 'level':level})
            except Exception:
                exc_type, exc_value, exc_tb = sys.exc_info()
                tbe = tb.TracebackException(
                    exc_type, exc_value, exc_tb,
                )
                log.error(logid+''.join(tbe.format()))

        pool.close()
        pool.join()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def judge_diff(raw, u, p, gs, ge, gstrand, ulim, cutoff, border, outdir, padding, queue=None, configurer=None, level=None):

    logid = scriptname+'.judge_diff: '
    try:
        if queue and level:
            configurer(queue, level)
        goi, chrom, strand, cons, reg, f, window, span = map(str,os.path.basename(raw).split(sep='_'))
        span = span.split(sep='.')[0]
        cs, ce = map(int, cons.split(sep='-'))
        ws, we = map(int, reg.split(sep='-'))

        cs = cs - ws #fit to window and make 0-based
        ce = ce - ws #fit to window and make 0-based closed

        if 0 > any([cs,ce,ws,we]):
            raise Exception('One of '+str([cs,ce,ws,we])+ ' lower than 0! this should not happen for '+','.join([goi, chrom, strand, cons, reg, f, window, span]))

        if gstrand != '-':
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
        log.debug(logid+'RAW: '+str(raw)+'\t'+str(noc))

        if abs(noc[ce]) > cutoff:
            uc = pl_to_array(u, ulim)  # This is the diffacc for unpaired constraint
            pc = pl_to_array(p, ulim)  # This is the diffacc for paired constraint

            log.debug(logid+'unpaired: '+str(u)+' and paired: '+str(p)+' Content: '+str(uc[ulim:ulim+10])+' test '+str(np.all(uc[ulim:ulim+10])))

            epsilon = 10 ** -50
            preaccu = noc+uc+epsilon
            preaccp = noc+pc+epsilon

            np.seterr(divide = 'ignore') #ignore 0 for LOGS
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

            np.seterr(divide = 'ignore') #ignore 0 for LOGS
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

            log.debug(logid+'WINDOWS: '+str.join(' ',map(str, [goi, strand, ws, cs, ce, we, str(cs+ws-1)+'-'+str(ce+ws), str(we-ce-1)+'-'+str(we-cs)])))

            for pos in range(len(noc)):
                if pos not in range(cs-padding+1-ulim, ce+padding+1+ulim):
                    if strand != '-':
                        gpos = pos + ws - ulim + 1  # already 0-based
                        gend = gpos + ulim  # 0-based half-open
                        gcst = cs+ws+1
                        gcen = ce+ws+2
                        gcons = str(gcst)+'-'+str(gcen)
                    else:
                        gpos = we - pos  # already 0-based
                        gend = gpos + ulim  # 0-based half-open
                        gcst = we-ce-1
                        gcen = we-cs
                        gcons = str(gcst)+'-'+str(gcen)

                    if border1 < uc[pos] and uc[pos] < border2:
                        if ce < pos:  # get distance up or downstream
                            dist = (pos - ce) * -1  # no -1 or we have 0 overlap
                        else:
                            dist = cs - pos

                        preacc = preaccu[pos] - epsilon
                        nrgdiff = nrgdiffu[pos]
                        kd = kdu[pos]
                        zscore = zscoresu[pos]

                        if not any([x is np.nan for x in [preacc,nrgdiff,kd,zscore]]):
                            out['u'].append('\t'.join([str(chrom), str(gpos), str(gend), str(goi) + '|' + str(cons) + '|' + str(gcons), str(uc[pos]), str(strand), str(dist), str(noc[pos]), str(preacc), str(nrgdiff), str(kd), str(zscore)]))

                    if border1 < pc[pos] and pc[pos] < border2:
                        if ce < pos:  # get distance up or downstream
                            dist = (pos - ce) * -1  # no -1 or we have 0 overlap
                        else:
                            dist = cs - pos

                        preacc = preaccu[pos] - epsilon
                        nrgdiff = nrgdiffu[pos]
                        kd = kdu[pos]
                        zscore = zscoresu[pos]

                        if not any([x is np.nan for x in [preacc,nrgdiff,kd,zscore]]):
                            out['p'].append('\t'.join([str(chrom), str(gpos), str(gend), str(goi) + '|' + str(cons) + '|' + str(gcons), str(pc[pos]), str(strand), str(dist), str(noc[pos]), str(preacc), str(nrgdiff), str(kd)]))

        savelists(out, outdir)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def savelists(out, outdir):

    logid = scriptname+'.savelist: '
    try:
        if len(out['u']) > 0:
            with gzip.open(os.path.abspath(os.path.join(outdir, 'Collection_unpaired.bed.gz')), 'ab') as o:
                o.write(bytes('\n'.join(out['u']),encoding='UTF-8'))
                o.write(bytes('\n',encoding='UTF-8'))
        if len(out['p']) > 0:
            with gzip.open(os.path.abspath(os.path.join(outdir, 'Collection_paired.bed.gz')), 'ab') as o:
                o.write(bytes('\n'.join(out['p']),encoding='UTF-8'))
                o.write(bytes('\n',encoding='UTF-8'))
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


def main(args):
    try:
        #  Logging configuration
        logdir = args.logdir
        ts = str(datetime.datetime.now().strftime("%Y%m%d_%H_%M_%S_%f"))
        logfile = str.join(os.sep,[os.path.abspath(logdir),scriptname+'_'+ts+'.log'])
        loglevel = args.loglevel

        makelogdir(logdir)
        makelogfile(logfile)

        queue = multiprocessing.Manager().Queue(-1)
        listener = multiprocessing.Process(target=listener_process, args=(queue, listener_configurer, logfile, loglevel))
        listener.start()

        worker_configurer(queue, loglevel)

        log.info(logid+'Running '+scriptname+' on '+str(args.procs)+' cores.')
        print('HERE')
        log.info(logid+'CLI: '+sys.argv[0]+' '+'{}'.format(' '.join( [shlex.quote(s) for s in sys.argv[1:]] )))

        screen_genes(queue, worker_configurer, loglevel, args.pattern, args.cutoff, args.border, args.ulimit, args.procs, args.roi, args.outdir, args.dir, args.genes, args.padding)

        queue.put_nowait(None)
        listener.join()


    except Exception:
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
        args=parseargs_collectpl()
        main(args)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

######################################################################
### CollectConsResults.py ends here
