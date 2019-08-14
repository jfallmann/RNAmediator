#!/usr/bin/env python3
### ConstraintPLFold.py ---
##
## Filename: ConstraintPLFold.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Sep  6 09:02:18 2018 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Wed Aug 14 14:24:12 2019 (+0200)
##           By: Joerg Fallmann
##     Update #: 151
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
############
#find ffmpeg executable
#import shutil
#plt.rcParams['animation.ffmpeg_path'] = shutil.which("ffmpeg")
#plt.rc('verbose', level='debug-annoying', fileo=sys.stdout)
#matplotlib.verbose.set_level("helpful")
#plt.rc('animation', html='html5')
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
streamlog = setup_multiprocess_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='WARNING')
log = setup_multiprocess_logger(name=scriptname, log_file='logs/'+scriptname, logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='INFO')

import argparse
import pprint
from io import StringIO
import time
import math
import gzip
import importlib
import multiprocessing
#from multiprocessing import Manager
import traceback as tb
from Randseq import createrandseq
#Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq
#numpy and matplolib and pyplot
import numpy as np
import matplotlib
from matplotlib import animation #, rc
import matplotlib.pyplot as plt
from random import choices, choice, shuffle # need this if tempprobing was choosen

def parseargs():
    parser = argparse.ArgumentParser(description='Calculate base pairing probs of given seqs or random seqs for given window size, span and region.')
    parser.add_argument("-s", "--sequence", type=str, help='Sequence to fold')
    parser.add_argument("-w", "--window", type=int, default=240, help='Size of window')
    parser.add_argument("-l", "--span", type=int, default=60, help='Length of bp span')
    parser.add_argument("-u", "--region", type=int, default=1, help='Length of region')
    parser.add_argument("-m", "--multi", type=int, default=4, help='Multiplyer for window expansion')
    parser.add_argument("-c", "--cutoff", type=float, default=-0.01, help='Only print prob greater cutoff')
    parser.add_argument("-r", "--unconstraint", type=str, default='STDOUT', help='Print output of unconstraint folding to file with this name')
    parser.add_argument("-n", "--unpaired", type=str, default='STDOUT', help='Print output of unpaired folding to file with this name')
    parser.add_argument("-p", "--paired", type=str, default='STDOUT', help='Print output of paired folding to file with this name')
    parser.add_argument("-e", "--length", type=int, default=100, help='Length of randseq')
    parser.add_argument("-g", "--gc", type=int, default=0, help='GC content, needs to be %%2==0 or will be rounded')
    parser.add_argument("-b", "--number", type=int, default=1, help='Number of random seqs to generate')
    parser.add_argument("-x", "--constrain", type=str, default='sliding', help='Region to constrain, either sliding window (default) or region to constrain (e.g. 1-10) or path to file containing regions following the naming pattern $fastaID_constraints e.g. Sequence1_constraints, if paired, the first entry of the file will become a fixed constraint and paired with all the others, choices = [off,sliding,temperature, tempprobe, the string file or a filename, paired, or simply 1-10,2-11 or 1-10;15-20,2-11;16-21 for paired or ono(oneonone),filename to use one line of constraint file for one sequence from fasta]')
    parser.add_argument("-y", "--conslength", type=int, default=0, help='Length of region to constrain for slidingwindow')
    parser.add_argument("-t", "--temprange", type=str, default='', help='Temperature range for structure prediction (e.g. 37-60)')
    parser.add_argument("-a", "--alphabet", type=str, default='AUCG', help='alphabet for random seqs')
    parser.add_argument("--plot", type=str, default='0', choices=['0','svg', 'png'], help='Create image of the (un-)constraint sequence, you can select the file format here (svg,png). These images can later on be animated with ImageMagick like `convert -delay 120 -loop 0 *.svg animated.gif`.')
    parser.add_argument("--save", type=int, default=1, help='Save the output as gz files')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory to write to')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processed to run this job with')
    parser.add_argument("--vrna", type=str, default='', help="Append path to vrna RNA module to sys.path")
    parser.add_argument("--pattern", type=str, default='', help="Helper var, only used if called from other prog where a pattern for files is defined")
    parser.add_argument("-v", "--verbosity", type=int, default=0, choices=[0, 1], help="increase output verbosity")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def preprocess(sequence, window, span, region, multi, unconstraint, unpaired, paired, length, gc, number, constrain, conslength, alphabet, plot, save, procs, vrna, temprange, outdir, verbosity=False, pattern=None, cutoff=None):

    logid = scriptname+'.preprocess: '

    #set path for output
    if outdir:
        if not os.path.isabs(outdir):
            outdir =  os.path.abspath(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    else:
        outdir = os.path.abspath(os.getcwd())

    if ( plot == '0' and not save):
        log.error(logid+'Neither plot nor save are active, this script will take a long time and produce nothing, please activate at least one of the two!')
        raise ValueError('Neither plot nor save are active, this script will take a long time and produce nothing, please activate at least one of the two!')

    mode = 'generic'
    if 'ono' == str(constrain.split(',')[0]):
        constrain = constrain.split(',')[1]
        mode = 'ono'
        constraintlist=[]
        if (os.path.isfile(constrain)):
            linewise = True
            if '.bed' in constrain:
                log.info(logid+'Parsing constraints from Bed '+constrain)
                if '.gz' in constrain:
                    f = gzip.open(constrain,'rt')
                else:
                    f = open(constrain,'rt')
                constraintlist = readConstraintsFromBed(f,linewise)
            elif '.csv' in constrain:
                if '.gz' in constrain:
                    f = gzip.open(constrain,'rt')
                else:
                    f = open(constrain,'rt')
                constraintlist = readConstraintsCSV(f,linewise)
            else:
                if '.gz' in constrain:
                    f = gzip.open(constrain,'rt')
                else:
                    f = open(constrain,'rt')
                constraintlist = readConstraintsFromGeneric(f,linewise)
            f.close()

        seq = parseseq(sequence)
        records = list(SeqIO.parse(seq, "fasta"))
        seqnr = 0
        # Create process pool with processes
        try:
            num_processes = procs or 1
            pool = multiprocessing.Pool(processes=num_processes-1, maxtasksperchild=1)
#           res = []
            for rec in records:
                #           with StringIO(records[seqnr].format("fasta")) as sseq:
                sseq = StringIO(records[seqnr].format("fasta"))
                constraint = constraintlist['lw'][seqnr]
                #           kwords = {'verbosity':verbosity, 'pattern':pattern, 'cutoff':cutoff, 'seqnr':seqnr, 'mode':mode, 'constraintlist' : 'Parsed'}
                try:
                    pool.apply_async(parafold, (sseq, window, span, region, multi, unconstraint, unpaired, paired, length, gc, number, constraint, conslength, alphabet, plot, save, procs, vrna, temprange, outdir, pattern, cutoff, seqnr))
                    seqnr += 1
                except Exception as err:
                    exc_type, exc_value, exc_tb = sys.exc_info()
                    tbe = tb.TracebackException(
                        exc_type, exc_value, exc_tb,
                    )
                    log.error(logid+''.join(tbe.format()))

            pool.close()
            pool.join()               #timeout
#           pool = None
#           res = []
        except Exception as err:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
                exc_type, exc_value, exc_tb,
                )
            log.error(logid+''.join(tbe.format()))

    else:
        try:
            fold(sequence, window, span, region, multi, unconstraint, unpaired, paired, length, gc, number, constrain, conslength, alphabet, plot, save, procs, vrna, temprange, outdir, verbosity=verbosity, pattern=pattern, cutoff=cutoff)
        except Exception as err:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
                exc_type, exc_value, exc_tb,
                )
            log.error(logid+''.join(tbe.format()))
    return 1

def fold(sequence, window, span, region, multi, unconstraint, unpaired, paired, length, gc, number, constrain, conslength, alphabet, plot, save, procs, vrna, temprange, outdir, verbosity=False, pattern=None, cutoff=None, seqnr=None, mode=None, constraintlist=None):

    logid = scriptname+'.fold: '
    #set path for VRNA lib if necessary
    if vrna:
        sys.path=[vrna] + sys.path
    try:
        seq = parseseq(sequence)

        global RNA
        RNA = importlib.import_module('RNA')
        globals().update(
        {n: getattr(RNA, n) for n in RNA.__all__}
            if hasattr(RNA, '__all__')
            else {k: v for (k, v) in RNA.__dict__.items() if not k.startswith('_')
            })
        md = RNA.md()
        md = None
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
        sys.exit()

    # Create process pool with processes
    num_processes = procs or 1
    pool = multiprocessing.Pool(processes=num_processes, maxtasksperchild=1)
#    logger = multiprocessing.log_to_stderr()   #LOGGING FOR DEBUG
#    logger.setLevel(logging.DEBUG)   #LOGGING FOR DEBUG

    for fa in SeqIO.parse(seq,'fasta'):
        fa.seq = str(fa.seq).upper()
        goi, chrom, strand = idfromfa(fa.id)

        if pattern and pattern not in goi:
            next
        else:
            log.info(logid+'Working on ' + goi + "\t" + fa.id)
##prepare plots
#       if plot != '0':
#           manager = Manager()
#           animations = manager.list()
            animations = []
            xvals = []
            for y in range(1,len(fa.seq)+1):
                xvals.append(y)
                xs = np.array(xvals)
                xvals = []
#set kT for nrg2prob and vice versa calcs
            kT = 0.61632077549999997
            #define data structures
            #           data = { 'bpp': [], 'up': [] }
            data = { 'up': [] }
            an = [np.nan]
            # We check if we need to fold the whole seq or just a region around the constraints
            if constrain == 'sliding' or constrain == 'temperature': # here we fold the whole seq
                #if not already available, run raw fold for whole sequence as we have a sliding window constraint
                check = str(goi)+'_'+str(chrom)+'_'+str(strand)+'_'+unconstraint+'_'+str(window)+'_'+str(span)+'.gz'
                if not ( os.path.isfile(check) ):
                    try:
                        res = [pool.apply_async(fold_unconstraint, args=(str(fa.seq), str(fa.id), region, window, span, unconstraint, save, outdir))]
                    except Exception as err:
                        exc_type, exc_value, exc_tb = sys.exc_info()
                        tbe = tb.TracebackException(
                            exc_type, exc_value, exc_tb,
                        )
                        log.error(logid+''.join(tbe.format()))
                        sys.exit()
                    data['up'] = res.get()
                else:
                    #The hard work of folding this has already been done, so we just read in the results
                    log.info(logid+'Found '+str(goi+'_'+unconstraint+'_'+str(window)+'_'+str(span)+'.gz')+', will read in data and not recalculate')
                    data['up'] = read_precalc_plfold(data['up'], check, str(fa.seq))

                an = up_to_array(data['up'],int(region),len(fa.seq))#convert to array for fast diff calc

            if constrain == 'temperature':
#                printlog('Calculating probs for temperature constraint '+temprange)
                ts, te = map(int,temprange.split('-'))

                try:
                    for temp in range(ts,te+1):
                    # Create the process, and connect it to the worker function
                        pool.apply_async(constrain_temp, args=(str(fa.id), str(fa.seq), temp, window, span, region, an, animations, xs, save, outdir, plot))
                except Exception as err:
                    exc_type, exc_value, exc_tb = sys.exc_info()
                    tbe = tb.TracebackException(
                        exc_type, exc_value, exc_tb,
                        )
                    log.error(logid+''.join(tbe.format()))

            elif constrain == 'tempprobe':
                probeargs = temprange.split(',')
                log.info(logid+'Calculating cutoff probs for temperature constraint '+probeargs[0])
                ts, te = map(int,probeargs[0].split('-'))
                if len(fa.seq) > window*multi:# If the sequence is very large and we only need a proxy of changes to get a border for the cutoff of CalcConsDiffs we randomly select 10 regions of size window of the sequence for folding
                    selection = ''
                    if probeargs[1]:# if we have constraints we choose windows around the constraints for tempprobing
                        for const in probeargs[1:]:
                            conss, conse = map(int,const.split('-'))
                            strt = conss-window*multi
                            endt = conse+window*multi+1
                            if strt < 0:
                                strt = 0
                            if endt > len(fa.seq):
                                endt = len(fa.seq)
                            selection = selection + str(fa.seq[strt:endt])
                    else:
                        for chosen in range(1,11): #we randomly choose 10 windows from the sequence
                            items = range(window,len(fa.seq)-window-1)
                            sel = choice(items)
                            selection = selection + (fa.seq[sel:sel+window])
                            fa.seq = Seq(str(selection))
                            #we check if we already have the raw seq folded
                check = str(goi)+'_'+str(chrom)+'_'+str(strand)+'_'+unconstraint+'_'+str(window)+'.gz'
                if ( os.path.isfile(check) ):
                    #The hard work of folding this has already been done, so we just read in the results
                    log.info(logid+'Found '+check+', will read in data and not recalculate')
                    data['up'] = read_precalc_plfold(data['up'], check, str(fa.seq))
                    #convert to array for fast diff calc
                    an = up_to_array(data['up'],int(region),len(fa.seq))

                if an[0] is np.nan or len(an) > len(fa.seq): #This means that the prob info was loaded from file or the raw sequence has not been folded yet, but we need a subsequence for the cutoff calculation, so we fold the subsequence
                    log.info(logid+'Recalculating at default temp with subseq')
                    try:
                        res = [pool.apply_async(fold_unconstraint, args=(str(fa.seq), str(fa.id), region, multi, window, span, unconstraint, save, outdir))]
                    except Exception as err:
                        exc_type, exc_value, exc_tb = sys.exc_info()
                        tbe = tb.TracebackException(
                            exc_type, exc_value, exc_tb,
                        )
                        log.error(logid+''.join(tbe.format()))

                    data['up'] = res.get()
                    an = up_to_array(data['up'],int(region),len(fa.seq))

                try:
                    for temp in range(ts,te+1):
# Create the process, and connect it to the worker function
                        pool.apply_async(constrain_temp, args=(str(fa.id), str(fa.seq), temp, window, span, region, multi, an, animations, xs, save, outdir, plot),error_callback=eprint)
                except Exception as err:
                    exc_type, exc_value, exc_tb = sys.exc_info()
                    tbe = tb.TracebackException(
                        exc_type, exc_value, exc_tb,
                        )
                    log.error(logid+''.join(tbe.format()))

            else:
                conslist = []
                if constraintlist is None or not constraintlist:
                    constraintlist=[]
                    if (os.path.isfile(os.path.abspath(constrain))):
                        constrain=os.path.abspath(constrain)
                        if '.bed' in constrain:
                            log.info(logid+'Parsing constraints for fold from Bed '+constrain)
                            if '.gz' in constrain:
                                f = gzip.open(constrain,'rt')
                            else:
                                f = open(constrain,'rt')
                            constraintlist = readConstraintsFromBed(f)
                        elif '.csv' in constrain:
                            if '.gz' in constrain:
                                f = gzip.open(constrain,'rt')
                            else:
                                f = open(constrain,'rt')
                            constraintlist = readConstraintsCSV(f)
                        else:
                            if '.gz' in constrain:
                                f = gzip.open(constrain,'rt')
                            else:
                                f = open(constrain,'rt')
                            constraintlist = readConstraintsFromGeneric(f)
                        f.close()

                if (os.path.isfile(os.path.abspath(constrain))):
                    try:
                        if goi in constraintlist:
                            conslist = constraintlist[goi]
                            log.info(logid+'Calculating probs for '+goi+' with constraint from file ' + constrain)
                        elif 'generic' in constraintlist:
                            conslist = constraintlist['generic']
                        elif 'NOCONS' in constraintlist:
                            conslist = constraintlist
                        elif constrain == 'sliding':
                            for start in range(1,len(fa.seq)-conslength+2):
                                end = start+conslength-1
                                conslist.append(str(start)+'-'+str(end))
                        elif ',' in constrain:
                            log.info(logid+'Calculating probs for constraint ' + constrain)
                            constraintlist = constrain.split(',')
                        else:
                            log.warning(logid+'Could not compute constraints for '+str(goi)+' May be missing in constraint file.')
                            continue
                    except Exception as err:
                        exc_type, exc_value, exc_tb = sys.exc_info()
                        tbe = tb.TracebackException(
                        exc_type, exc_value, exc_tb,
                        )
                        log.error(logid+''.join(tbe.format()))

                elif constrain == 'file' or constrain == 'paired':
                    log.info(logid+'Calculating probs for constraint from file ' + str(goi + '_constraints'))
                    with open(goi+'_constraints','rt') as o:
                        for line in o:
                            conslist.append(line.rstrip())
                    o.close()
                elif constrain == 'none':
                    conslist = ['NOCONS']
                elif constrain == 'sliding':
                    for start in range(1,len(fa.seq)-conslength+2):
                        end = start+conslength-1
                        conslist.append(str(start)+'-'+str(end))
                else:
                    log.info(logid+'Calculating probs for constraint ' + constrain)
                    conslist = constrain.split(',')

                log.debug(conslist)

                for entry in conslist:
                    if entry == 'NOCONS': # in case we just want to fold the sequence without constraints at all
                        try:
                            res = [pool.apply_async(fold_unconstraint, args=(str(fa.seq), str(fa.id), region, window, span, multi, unconstraint, save, outdir))]
                        except Exception as err:
                            exc_type, exc_value, exc_tb = sys.exc_info()
                            tbe = tb.TracebackException(
                                exc_type, exc_value, exc_tb,
                            )
                            log.error(logid+''.join(tbe.format()))

                        data['up'] = res.get()

                    else:
                        # we now have a list of constraints and for the raw seq comparison we only need to fold windows around these constraints
                        # In case we want to constrain pairwise
                        fstart, fend = [None,None]
                        if constrain == 'paired' or ':' in entry:
                            [fstart, fend], [start, end] = [[int(x) for x in cn.split('-',1)] for cn in entry.split(':',1)]
                            cons = str(fstart)+'-'+str(fend)+':'+str(start)+'-'+str(end)
                            if start < 0 or fstart < 0 or end > len(fa.seq) or fend > len(fa.seq):
                                log.warning(logid+'Constraint out of sequence bounds! skipping! '+','.join([len(fa.seq),str(start)+'-'+str(end),str(fstart)+'-'+str(fend)]) )
                                continue

                        else:
                            start,end = map(int,entry.split('-',1))
                            tostart, toend = expand_window(start, end, window, multi, len(fa.seq))
                            cons = str(start)+'-'+str(end)+'_'+str(tostart)+'-'+str(toend)

                            if start < 0 or end > len(fa.seq):
                                log.warning(logid+'Constraint out of sequence bounds! skipping! '+','.join(map(str,[len(fa.seq),str(start)+'-'+str(end)])) )
                                continue

                        if checkexisting(str(fa.id), paired, unpaired, cons, region, window, span, outdir):
                            log.warning(logid+str(cons)+' Exists for '+str(fa.id)+'! Skipping!')
                            continue

                        log.info(logid+'Calculating constraint\t' + entry)
                        const = np.array([start, end])

                        data = { 'up': [] }
                        an = None

                        if fstart is not None and fend is not None:
                            log.info(logid+'Constraining to '+str(fstart) + ' and ' + str(fend))
                            goi, chrom, strand = idfromfa(fa.id)

                            try:
                                pool.apply_async(constrain_seq_paired, args=(str(fa.id), str(fa.seq), fstart, fend, start, end, conslength, const, cons, window, span, region, multi, animations, xs, paired, unpaired, save, outdir, plot, data, an, unconstraint))
                            except Exception as err:
                                exc_type, exc_value, exc_tb = sys.exc_info()
                                tbe = tb.TracebackException(
                                    exc_type, exc_value, exc_tb,
                                )
                                log.error(logid+''.join(tbe.format()))
                        else:
                            try:
                                pool.apply_async(constrain_seq, args=(str(fa.id), str(fa.seq), start, end, conslength, const, cons, window, span, region, multi, animations, xs, paired, unpaired, save, outdir, plot, data, an, unconstraint))
                            except Exception as err:
                                exc_type, exc_value, exc_tb = sys.exc_info()
                                tbe = tb.TracebackException(
                                    exc_type, exc_value, exc_tb,
                                )
                                log.error(logid+''.join(tbe.format()))

#       results = [p.get for p in res]   #collecting results, shouls not return values but make processes wait for workers to finish, Slows down everything, hopefully removable
    try:
        pool.close()
        pool.join()       #timeout
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

    log.info(logid+"DONE: output in: " + str(outdir))
    return 1

def parafold(sequence, window, span, region, multi, unconstraint, unpaired, paired, length, gc, number, constrain, conslength, alphabet, plot, save, procs, vrna, temprange, outdir, pattern, cutoff, seqnr):

    logid = scriptname+'.parafold: '
    seq = sequence
    #set path for VRNA lib
    if vrna:
        sys.path=[vrna] + sys.path
    else:
        sys.path=["/scratch/fall/VRNA/243_alpha2/lib/python3.6/site-packages"] + sys.path
    try:
        global RNA
        RNA = importlib.import_module('RNA')
        globals().update(
            {n: getattr(RNA, n) for n in RNA.__all__}
            if hasattr(RNA, '__all__')
            else {k: v for (k, v) in RNA.__dict__.items() if not k.startswith('_')
        })
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

    try:
        for fa in SeqIO.parse(seq,'fasta'):
            fa.seq = str(fa.seq).upper()
            goi, chrom, strand = idfromfa(fa.id)

            if pattern and pattern not in goi:
                next
            else:
                log.info(logid+'Working on ' + goi)

            # if plot != '0':
            # manager = Manager() ###AssertionError: daemonic processes are not allowed to have children
            # animations = manager.list()
            animations = ()     # we have to add this for now
            xvals = []
            for y in range(1,len(fa.seq)+1):
                xvals.append(y)
            xs = np.array(xvals)
            xvals = []
            kT = 0.61632077549999997
            #define data structures
            #           data = { 'bpp': [], 'up': [] }
            data = { 'up': [] }
            an = [np.nan]
            num_processes = procs or 1
            conslist = []
            conslist.append(constrain)

            for entry in conslist:
                if entry == 'NOCONS': # in case we just want to fold the sequence without constraints at all
                    try:
                        md = RNA.md()
                        md.max_bp_span = span
                        md.window_size = window
                        # create new fold_compound object
                        fc = RNA.fold_compound(str(seq), md, RNA.OPTION_WINDOW)
                        # call prop window calculation
                        fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data)
                        #fc.probs_window(region, RNA.PROBS_WINDOW_BPP, bpp_callback, data)
                    except Exception as err:
                        exc_type, exc_value, exc_tb = sys.exc_info()
                        tbe = tb.TracebackException(
                            exc_type, exc_value, exc_tb,
                            )
                        log.error(logid+''.join(tbe.format()))

                    if save:
                        write_unconstraint(str(fa.id), str(fa.seq), unconstraint, data, int(region), str(window), str(span), outdir)

                else:
                    # we now have a list of constraints and for the raw seq comparison we only need to fold windows around these constraints
                    # In case we want to constrain pairwise
                    data = { 'up': [] }
                    start,end = map(int,entry.split('-',1))

                    if start < 0 or end > len(fa.seq):
                        log.warning(logid+'Constraint out of sequence bounds! skipping! '+','.join([len(fa.seq),str(start)+'-'+str(end)]))
                        continue

                    log.info(logid+'Calculating constraint\t' + entry)
                    tostart, toend = expand_window(start, end, window, multi, len(fa.seq))

                    seqtofold = str(fa.seq[tostart:toend]).upper()
                    const = np.array([start, end])

                    log.info(logid+'Constraint: ' + str(entry) + '\tSeqlength: ' + str(len(fa.seq)) + '\tFoldlength: ' + str(len(seqtofold)) + ' : ' + str(tostart) + " - " + str(toend))
                    md = RNA.md()
                    md.max_bp_span = span
                    md.window_size = window

                    fc = RNA.fold_compound(str(seqtofold), md, RNA.OPTION_WINDOW)
                    data_t = { 'up': [] }
                    fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_t)
                    data['up'].extend(data_t['up'][:])
                    data_t = None # remove unneeded list
                    an = up_to_array(data['up'],int(region),len(seqtofold))

                    try:
                        constrain_seq(str(fa.id), str(fa.seq), start, end, conslength, const, cons, window, span, region, multi, animations, xs, paired, unpaired, save, outdir, plot, data, an)
                    except Exception as err:
                        exc_type, exc_value, exc_tb = sys.exc_info()
                        tbe = tb.TracebackException(
                            exc_type, exc_value, exc_tb,
                        )
                        log.error(logid+''.join(tbe.format()))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

    log.info(logid+"DONE: output in: " + str(outdir))
    return 1

##### Functions #####
def constrain_seq(sid, seq, start, end, conslength, const, cons, window, span, region, multi, animations, xs, paired, unpaired, save, outdir, plot, data, an=None, unconstraint=None):
    #   DEBUGGING
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    logid = scriptname+'.constrain_seq: '
    try:
        goi, chrom, strand = idfromfa(sid)

        #for all constraints we now extract subsequences to compare against
        #we no longer fold the whole raw sequence but only the constraint region +- window size
        tostart, toend = expand_window(start, end, window, multi, len(seq))
        seqtofold = str(seq[tostart:toend]) ###TEST

        cons = str('-'.join([str(start),str(end)])+'_'+'-'.join([str(tostart),str(toend)]))

        if start < 0 or end > len(seq):
            log.warning(logid+'Constraint out of sequence bounds! skipping! '+','.join([len(seq),str(start)+'-'+str(end)]) )
            return

        if checkexisting(sid, paired, unpaired, cons, region, window, span, outdir):
            log.warning(logid+str(cons)+' Existst for '+str(sid)+'! Skipping!')
            return

        #refresh model details
        md = RNA.md()
        md.max_bp_span = span
        md.window_size = window

        # create new fold_compound objects
        fc_p = RNA.fold_compound(seqtofold, md, RNA.OPTION_WINDOW)
        fc_u = RNA.fold_compound(seqtofold, md, RNA.OPTION_WINDOW)

        log.debug(''.join(map(str,[logid, sid, region, seqtofold, cons, start, end, tostart, start-tostart, end-tostart+1])))

        #enforce paired
        for x in range(start-tostart, end-tostart+1):
            fc_p.hc_add_bp_nonspecific(x,0) #0 means without direction  ( $ d < 0 $: pairs upstream, $ d > 0 $: pairs downstream, $ d == 0 $: no direction)

        #enforce unpaired
        for x in range(start-tostart, end-tostart+1):
            fc_u.hc_add_up(x)

        #new data struct
        data_pn = {'up': []}
        data_un = {'up': []}

        fc_p.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_pn)
        fc_u.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_un)

        #we now fill the list with values from the unconstraint sequence to get the same length
        data_p = {'up': []}
        data_u = {'up': []}

        data_p['up'].extend(data_pn['up'][:])
        data_pn = None # remove unneeded list
        data_u['up'].extend(data_un['up'][:])
        data_un = None # remove unneeded list

        au = up_to_array(data_u['up'])#,region,len(seqtofold))
        ap = up_to_array(data_p['up'])#,region,len(seqtofold))

        if not an or len(an) < 1 or len(data['up']) < 1:
            data['up'] = fold_unconstraint(str(seqtofold), sid, region, window, span, unconstraint, save, outdir, cons)
            an = up_to_array(data['up'])#,int(region),len(seqtofold))

        if not np.array_equal(an, au):
            diff_nu = an - au
        else:
            log.info(logid+'No influence on Structure with unpaired constraint at ' + cons)
            diff_nu = None
        if not np.array_equal(an, ap):
            diff_np = an - ap
        else:
            log.info(logid+'No influence on Structure with paired constraint at ' + cons)
            diff_np = None

        if plot == 'svg' or plot == 'png':
            plot_data(str(sid), seqtofold, an, au, ap, const, xs, cons, plot, outdir)

        if save:
            write_constraint(str(sid), seqtofold, paired, unpaired, data_u, data_p, cons, int(region), diff_nu, diff_np, str(window), str(span), outdir)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

    return 1

def constrain_seq_paired(sid, seq, fstart, fend, start, end, conslength, const, cons, window, span, region, multi, animations, xs, paired, unpaired, save, outdir, plot, data, an, unconstraint):
    #   DEBUGGING
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    logid = scriptname+'.constrain_seq_paired: '
    try:
        #we no longer fold the whole sequence but only the constraint region +- window size
        tostart, toend = expand_window(start, end, window, multi, len(seq))
        seqtofold = str(seq[tostart:toend]) ###TEST

        if start < 0 or end > len(seq) or fstart < 0 or fend > len(seq):
            log.warning(logid+'Constraint out of sequence bounds! skipping! '+','.join([len(seq), str(start)+'-'+str(end), str(fstart)+'-'+str(fend)]) )
            return

        if checkexisting(sid, paired, unpaired, cons, region, window, span, outdir):
            log.warning(logid+str(cons)+' Existst for '+str(sid)+'! Skipping!')
            return

        cons = str('-'.join([str(start),str(end)])+'_'+'-'.join([str(tostart),str(toend)]))

        #refresh model details
        md = RNA.md()
        md.max_bp_span = span
        md.window_size = window

        # create new fold_compound objects
        fc_p = RNA.fold_compound(seqtofold, md, RNA.OPTION_WINDOW)
        fc_u = RNA.fold_compound(seqtofold, md, RNA.OPTION_WINDOW)

        #enforce paired
        for x in range(start-tostart, end-tostart+1):
            fc_p.hc_add_bp_nonspecific(x,0) #0 means without direction  ( $ d < 0 $: pairs upstream, $ d > 0 $: pairs downstream, $ d == 0 $: no direction)
        #enforce paired
        for x in range(fstart-tostart, fend-tostart+1):
            fc_p.hc_add_bp_nonspecific(x,0) #0 means without direction  ( $ d < 0 $: pairs upstream, $ d >

        #enforce unpaired
        for x in range(start-tostart, end-tostart+1):
            fc_u.hc_add_up(x)
        #enforce unpaired
        for x in range(fstart-tostart, fend-tostart+1):
            fc_u.hc_add_up(x)

        #new data struct
        data_pn = {'up': []}
        data_un = {'up': []}

        fc_p.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_pn)
        fc_u.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_un)

        #we now fill the list with values from the unconstraint sequence to get the same length
        data_p = {'up': []}
        data_u = {'up': []}

        data_p['up'].extend(data_pn['up'][:])
        data_pn = None # remove unneeded list
        data_u['up'].extend(data_un['up'][:])
        data_un = None # remove unneeded list

        au = up_to_array(data_u['up'],region,len(seqtofold))
        ap = up_to_array(data_p['up'],region,len(seqtofold))

        if not an or len(an) < 1 or len(data['up']) < 1:
            data['up'] = fold_unconstraint(str(seqtofold), sid, region, window, span, unconstraint, save, outdir, cons)
            an = up_to_array(data['up'],int(region),len(seqtofold))

        if not np.array_equal(an, au):
            diff_nu = an - au
        else:
            log.info(logid+'No influence on Structure with unpaired constraint at ' + cons)
            diff_nu = None
        if not np.array_equal(an, ap):
            diff_np = an - ap
        else:
            log.info(logid+'No influence on Structure with paired constraint at ' + cons)
            diff_np = None

        if plot == 'svg' or plot == 'png':
            plot_data(sid, seqtofold, an, au, ap, const, xs, cons, plot, outdir)

        if save:
            write_constraint(sid, seqtofold, paired, unpaired, data_u, data_p, cons, int(region), diff_nu, diff_np, str(window), str(span), outdir)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

    return 1

def constrain_temp(sid, seq, temp, window, span, region, multi, an, animations, xs, save, outdir, plot):

    #refresh model details
    md = RNA.md()
    md.max_bp_span = span
    md.window_size = window
    #set temperature
    md.temperature = temp
    #create new fold_compound objects
    fc_t = RNA.fold_compound(str(seq), md, RNA.OPTION_WINDOW)
    #new data struct
    data_t = {'up': []}

    fc_t.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_t)

    at = up_to_array(data_t['up'],region,len(seq))

    diff_nt = an - at

#####Set temp, and rewrite save and plot for temp
    if save:
        write_temp(sid, seq, str(temp), data_t, int(region), diff_nt, str(window), str(span), outdir)

    if plot == 'svg' or plot == 'png':
        plot_temp(sid, seq, at, temp, xs, plot, outdir)

    return 1

def write_unconstraint(sid, seq, unconstraint, data, region, window, span, outdir, rawentry=None):

    try:
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir,goi)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

    try:
        gr = str(sid.split(':')[3].split('(')[0])
    except:
        gr = 'na'

    try:
        if unconstraint != 'STDOUT':
            if not os.path.exists(temp_outdir):
                os.makedirs(temp_outdir)
            if rawentry:
                if not os.path.exists(os.path.join(temp_outdir,str(goi+'_'+chrom+'_'+strand+'_'+rawentry+'_'+unconstraint+'_'+window+'_'+str(span)+'.gz'))):
                    with gzip.open(os.path.join(temp_outdir,goi+'_'+chrom+'_'+strand+'_'+rawentry+'_'+unconstraint+'_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                        out = print_up(data['up'],len(seq),region)
                        if out and len(out)>1:
                            o.write(bytes(out,encoding='UTF-8'))
                        else:
                            print("No output produced "+sid)
            else:
                if not os.path.exists(os.path.join(temp_outdir,str(goi+'_'+chrom+'_'+strand+'_'+str(gr)+'_'+unconstraint+'_'+window+'_'+str(span)+'.gz'))):
                    with gzip.open(os.path.join(temp_outdir,goi+'_'+chrom+'_'+strand+'_'+str(gr)+'_'+unconstraint+'_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                        out = print_up(data['up'],len(seq),region)
                        if out and len(out)>1:
                            o.write(bytes(out,encoding='UTF-8'))
        else:
            print (print_up(data['up'],len(seq),region))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
    return 1

def write_constraint(sid, seq, paired, unpaired, data_u, data_p, constrain, region, diff_nu, diff_np, window, span, outdir):

    try:
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir,goi)
#print outputs to file or STDERR
        if paired != 'STDOUT':
            if not os.path.exists(temp_outdir):
                os.makedirs(temp_outdir)
            if not os.path.exists(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_'+paired+'_'+window+'_'+str(span)+'.gz')):
                with gzip.open(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_'+paired+'_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                    out = print_up(data_p['up'],len(seq),region)
                    if out and len(out)>1:
                        o.write(bytes(out,encoding='UTF-8'))
                    else:
                        eprint("No output produced "+sid)
        else:
            print(print_up(data_p['up'],len(seq),region))

        if unpaired != 'STDOUT':
            if not os.path.exists(temp_outdir):
                os.makedirs(temp_outdir)
            if not os.path.exists(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_'+unpaired+'_'+window+'_'+str(span)+'.gz')):
                with gzip.open(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_'+unpaired+'_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                    out = print_up(data_u['up'],len(seq),region)
                    if out and len(out)>1:
                        o.write(bytes(out,encoding='UTF-8'))
                    else:
                        eprint("No output produced "+sid)
        else:
            print(print_up(data_u['up'],len(seq),region))

        if diff_nu.any():
            if unpaired != 'STDOUT':
                if not os.path.exists(temp_outdir):
                    os.makedirs(temp_outdir)
                if not os.path.exists(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_diffnu_'+window+'_'+str(span)+'.gz')):
                    with gzip.open(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_diffnu_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                        printdiff(diff_nu,o)
            else:
                npprint(diff_nu)

        if diff_np.any():
            if unpaired != 'STDOUT':
                if not os.path.exists(temp_outdir):
                    os.makedirs(temp_outdir)
                if not os.path.exists(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_diffnp_'+window+'_'+str(span)+'.gz')):
                    with gzip.open(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_diffnp_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                        printdiff(diff_np,o)
            else:
                npprint(diff_np)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
    return 1

def prepare_write_ucons(sid, seq, unconstraint, data, region, window, span, outdir, rawentry=None):

    try:
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir,goi)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

    try:
        gr = str(sid.split(':')[3].split('(')[0])
    except:
        gr = 'na'

    try:
        if unconstraint != 'STDOUT':
            if not os.path.exists(temp_outdir):
                os.makedirs(temp_outdir)
            if rawentry:
                if not os.path.exists(os.path.join(temp_outdir,str(goi+'_'+chrom+'_'+strand+'_'+unconstraint+'_'+rawentry+'_'+window+'_'+str(span)+'.gz'))):
                    with gzip.open(os.path.join(temp_outdir,goi+'_'+chrom+'_'+strand+'_'+unconstraint+'_'+rawentry+'_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                        out = print_up(data['up'],len(seq),region)
                        if out and len(out)>1:
                            o.write(bytes(out,encoding='UTF-8'))
                        else:
                            log.error(logid+"No output produced "+sid)
            else:
                if not os.path.exists(os.path.join(temp_outdir,str(goi+'_'+chrom+'_'+strand+'_'+unconstraint+'_'+str(gr)+'_'+window+'_'+str(span)+'.gz'))):
                    with gzip.open(os.path.join(temp_outdir,goi+'_'+chrom+'_'+strand+'_'+unconstraint+'_'+str(gr)+'_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                        out = print_up(data['up'],len(seq),region)
                        if out and len(out)>1:
                            o.write(bytes(out,encoding='UTF-8'))
        else:
            print (print_up(data,len(seq),region))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
    return 1

def prepare_write_cons(sid, seq, paired, unpaired, data_u, data_p, constrain, region, diff_nu, diff_np, window, span, outdir):

    try:
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir,goi)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
#print outputs to file or STDERR
    try:
        if paired != 'STDOUT':
            if not os.path.exists(temp_outdir):
                os.makedirs(temp_outdir)
            if not os.path.exists(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_'+paired+'_'+window+'_'+str(span)+'.gz')):
                with gzip.open(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_'+paired+'_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                    out = print_up(data_p['up'],len(seq),region)
                    if out and len(out)>1:
                        o.write(bytes(out,encoding='UTF-8'))
                    else:
                        log.error(logid+"No output produced "+sid)
        else:
            print(print_up(data_p['up'],len(seq),region))

        if unpaired != 'STDOUT':
            if not os.path.exists(temp_outdir):
                os.makedirs(temp_outdir)
            if not os.path.exists(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_'+unpaired+'_'+window+'_'+str(span)+'.gz')):
                with gzip.open(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_'+unpaired+'_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                    out = print_up(data_u['up'],len(seq),region)
                    if out and len(out)>1:
                        o.write(bytes(out,encoding='UTF-8'))
                    else:
                        eprint("No output produced "+sid)
        else:
            print(print_up(data_u['up'],len(seq),region))

        if diff_nu.any():
            if unpaired != 'STDOUT':
                if not os.path.exists(temp_outdir):
                    os.makedirs(temp_outdir)
                if not os.path.exists(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_diffnu_'+window+'_'+str(span)+'.gz')):
                    with gzip.open(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_diffnu_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                        printdiff(diff_nu,o)
            else:
                npprint(diff_nu)

        if diff_np.any():
            if unpaired != 'STDOUT':
                if not os.path.exists(temp_outdir):
                    os.makedirs(temp_outdir)
                if not os.path.exists(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_diffnp_'+window+'_'+str(span)+'.gz')):
                    with gzip.open(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+constrain+'_diffnp_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                        printdiff(diff_np,o)
            else:
                npprint(diff_np)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
    return 1


def write_temp(sid, seq, temp, data, region, diff, window, span, outdir):

    logid = scriptname+'.write_temp: '
    #print outputs to file or STDERR
    try:
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir,goi)
        if not os.path.exists(temp_outdir):
            os.makedirs(temp_outdir)
        if not os.path.exists(os.path.join(temp_outdir,'TempCons_'+goi+'_'+chrom+'_'+strand+'_'+temp+'_temp_'+window+'_'+str(span)+'.gz')):
            with gzip.open(os.path.join(temp_outdir,'TempCons_'+goi+'_'+chrom+'_'+strand+'_'+temp+'_temp_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                o.write(bytes(print_up(data['up'],len(seq),region),encoding='UTF-8'))
            if not os.path.exists(os.path.join(temp_outdir,'TempCons_'+goi+'_'+chrom+'_'+strand+'_'+temp+'_difft_'+window+'_'+str(span)+'.gz')):
                with gzip.open(os.path.join(temp_outdir,'TempCons_'+goi+'_'+chrom+'_'+strand+'_'+temp+'_difft_'+window+'_'+str(span)+'.gz'), 'wb') as o:
                    npprint(diff,o)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
    return 1

def checkexisting(sid, paired, unpaired, cons, region, window, span, outdir):

    logid = scriptname+'.checkexisting: '
    try:
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir,goi)

        if os.path.exists(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+cons+'_'+paired+'_'+str(window)+'_'+str(span)+'.gz')) and os.path.exists(os.path.join(temp_outdir,'StruCons_'+goi+'_'+chrom+'_'+strand+'_'+cons+'_'+unpaired+'_'+str(window)+'_'+str(span)+'.gz')):
            return True
        else:
            return False
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
    return 1

def bpp_callback(v, v_size, i, maxsize, what, data):

    logid = scriptname+'.bpp_callback: '
    try:
        if what & RNA.PROBS_WINDOW_BPP:
            data['bpp'].extend([{'i': i, 'j': j, 'p': p} for j, p in enumerate(v) if (p is not None)])# and (p >= 0.01)])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def up_callback(v, v_size, i, maxsize, what, data):

    logid = scriptname+'.up_callback: '
    try:
        if what & RNA.PROBS_WINDOW_UP:
            data['up'].extend([v])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def fold_unconstraint(seq, id, region, window, span, unconstraint, save, outdir, rawentry=None):

    logid = scriptname+'.fold_unconstraint: '
    data = { 'up': [] }
    try:
        md = RNA.md()
        md.max_bp_span = span
        md.window_size = window
        # create new fold_compound object
        fc = RNA.fold_compound(str(seq), md, RNA.OPTION_WINDOW)
        # call prop window calculation
        fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data)
        #fc.probs_window(region, RNA.PROBS_WINDOW_BPP, bpp_callback, data)
        if save:
            write_unconstraint(str(id), str(seq), unconstraint, data, int(region), str(window), str(span), outdir, rawentry)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


    return data['up']

def expand_window(start, end, window, multiplyer, seqlen):

    logid = scriptname+'.expand_window: '
    try:
        tostart = start - multiplyer*window
        if tostart < 0:
            tostart = 0
        toend = end + multiplyer*window + 1
        if toend > seqlen:
            toend = seqlen
        return [tostart, toend]
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
        preprocess(args.sequence, args.window, args.span, args.region, args.multi, args.unconstraint, args.unpaired, args.paired, args.length, args.gc, args.number, args.constrain, args.conslength, args.alphabet, args.plot, args.save, args.procs, args.vrna, args.temprange, args.outdir, args.verbosity, args.pattern, args.cutoff)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

    # ConstraintPLFold.py ends here
