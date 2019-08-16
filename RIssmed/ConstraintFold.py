#!/usr/bin/env python
### ConstraintFold.py ---
##
## Filename: ConstraintFold.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Sep  6 09:02:18 2018 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Fri Aug 16 09:36:56 2019 (+0200)
##           By: Joerg Fallmann
##     Update #: 344
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
global streamlog, log           # global to ensure that later manipulation of loglevel is possible
streamlog = setup_multiprocess_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='WARNING')
log = setup_multiprocess_logger(name=scriptname, log_file='logs/'+scriptname, logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='WARNING')

#other modules
import argparse
import pprint
from io import StringIO
import time
import math
import gzip
import copy
import importlib
import multiprocessing
from multiprocessing import Manager
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
import collections

def parseargs():
    parser = argparse.ArgumentParser(description='Calculate base pairing probs of given seqs or random seqs for given window size, span and region.')
    parser.add_argument("-s", "--sequence", type=str, help='Sequence to fold')
    parser.add_argument("-w", "--window", type=int, default=None, help='Size of window')
    parser.add_argument("-l", "--span", type=int, default=None, help='Maximum bp span')
    parser.add_argument("-c", "--cutoff", type=float, default=-0.01, help='Only print prob greater cutoff')
    parser.add_argument("-r", "--unconstraint", type=str, default='STDOUT', help='Print output of unconstraint folding to file with this name, use "STDOUT" to print to STDOUT (default)')# or "add" to append to paired/unpaired output')
    parser.add_argument("-n", "--unpaired", type=str, default='STDOUT', help='Print output of unpaired folding to file with this name')
    parser.add_argument("-p", "--paired", type=str, default='STDOUT', help='Print output of paired folding to file with this name')
    parser.add_argument("-e", "--length", type=int, default=100, help='Length of randseq')
    parser.add_argument("-g", "--gc", type=int, default=0, help='GC content, needs to be %%2==0 or will be rounded')
    parser.add_argument("-b", "--number", type=int, default=1, help='Number of random seqs to generate')
    parser.add_argument("-x", "--constrain", type=str, default='sliding', help='Region to constrain, either sliding window (default) or region to constrain (e.g. 1-10) or path to file containing regions following the naming pattern $fastaID_constraints, if paired, the first entry of the file will become a fixed constraint and paired with all the others, e.g. Sequence1_constraints, choices = [off,sliding,temperature, tempprobe, file, paired, or simply 1-10,2-11 or 1-10;15-20,2-11:16-21 for paired]')
    parser.add_argument("-y", "--conslength", type=int, default=0, help='Length of region to constrain for slidingwindow')
    parser.add_argument("-t", "--temprange", type=str, default='', help='Temperature range for structure prediction (e.g. 37-60)')
    parser.add_argument("-a", "--alphabet", type=str, default='AUCG', help='alphabet for random seqs')
    parser.add_argument("--plot", type=str, default='0', choices=['0','svg', 'png'], help='Create image of the (un-)constraint sequence, you can select the file format here (svg,png). These images can later on be animated with ImageMagick like `convert -delay 120 -loop 0 *.svg animated.gif`.')
    parser.add_argument("--save", type=int, default=1, help='Save the output as gz files')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory to write to')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processed to run this job with')
    parser.add_argument("--vrna", type=str, default='', help="Append path to vrna RNA module to sys.path")
    parser.add_argument("--pattern", type=str, default='', help="Helper var, only used if called from other prog where a pattern for files is defined")
    parser.add_argument("-i", "--genes", type=str, help='Genomic coordinates bed for genes, either standard bed format or AnnotateBed.pl format')
    parser.add_argument("-v", "--verbosity", type=int, default=0, choices=[0, 1], help="increase output verbosity")
    parser.add_argument("--loglevel", type=str, default='WARNING', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def preprocess(sequence, window, span, unconstraint, unpaired, paired, length, gc, number, constrain, conslength, alphabet, plot, save, procs, vrna, temprange, outdir, genes, verbosity=False, pattern=None, cutoff=None):

    logid = scriptname+'.preprocess: '
    try:
        #set path for output
        if outdir:
            if not os.path.isabs(outdir):
                outdir =  os.path.abspath(outdir)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
        else:
            outdir = os.path.abspath(os.getcwd())

        if ( plot == '0' and not save):
            raise ValueError('Neither plot nor save are active, this script will take a long time and produce nothing, please activate at least one of the two!')
        genecoords = parse_annotation_bed(genes) #get genomic coords to print to bed later, should always be just one set of coords per gene

        fold(sequence, window, span, unconstraint, unpaired, paired, length, gc, number, constrain, conslength, alphabet, plot, save, procs, vrna, temprange, outdir, genecoords, verbosity=verbosity, pattern=pattern, cutoff=cutoff)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def fold(sequence, window, span, unconstraint, unpaired, paired, length, gc, number, constrain, conslength, alphabet, plot, save, procs, vrna, temprange, outdir, genecoords, verbosity=False, pattern=None, cutoff=None):

    logid = scriptname+'.fold: '
    try:
        seq = parseseq(sequence)

        #set path for VRNA lib if necessary
        if vrna:
            sys.path=[vrna] + sys.path
        try:
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
            return

        # Create process pool with processes
        num_processes = procs or 1
        pool = multiprocessing.Pool(processes=num_processes, maxtasksperchild=1)
    #    logger = multiprocessing.log_to_stderr()   #LOGGING FOR DEBUG
    #    logger.setLevel(logging.DEBUG)

        if ( plot == '0' and not save):
            raise ValueError('Neither plot nor save are active, this script will take a long time and produce nothing, please activate at least one of the two!')

        ##prepare plots
        #       if plot != '0':
        manager = Manager()
        animations = manager.list()
        outdict = collections.defaultdict()  # Does not throw keyerror on missing key

        # constraints
        if constrain == 'temperature':
            log.info(logid+'Calculating probs for temperature constraint '+temprange)
            # Create process pool with processes
            num_processes = procs
            pool = multiprocessing.Pool(processes=num_processes)
            processes = []

            ts, te = map(int,temprange.split('-'))

            for fa in SeqIO.parse(seq,'fasta'):
                goi, chrom, strand = idfromfa(fa.id)

                if pattern and pattern not in goi:
                    next
                else:
                    log.info(logid+'Working on ' + goi)
                    xvals = []
                    for y in range(1,len(fa.seq)+1):
                        xvals.append(y)
                    xs = np.array(xvals)
                    xvals = []

                for temp in range(ts,te+1):
                    # Create the process, and connect it to the worker function
                    pool.apply_async(constrain_temp, args=(fa, temp, window, span, an, animations, xs, save, outdir, plot))

            pool.close()
            pool.join()

        else:
            #Create process pool with processes
            num_processes = procs
            pool = multiprocessing.Pool(processes=num_processes)
            conslist = list()
            constraintlist= list()

            if (os.path.isfile(constrain)):
                if '.bed' in constrain:
                    log.info(logid+'Parsing constraints for fold from Bed '+constrain)
                    if '.gz' in constrain:
                        f = gzip.open(constrain,'rt')
                    else:
                        f = open(constrain,'rt')
                    if any(x in constrain for x in ['paired','Paired']):
                        log.info(logid+'Reading paired constraints')
                        constraintlist = readPairedConstraintsFromBed(f)
                    else:
                        log.info(logid+'Reading single constraints')
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
            elif constrain == 'file' or constrain == 'paired':
                log.info(logid+'Calculating probs for constraint from file ' + str(goi + '_constraints'))
                with open(goi+'_constraints','rt') as o:
                    for line in o:
                        conslist.append(line.rstrip())
            elif constrain == 'none':
                constraintlist = ['NOCONS']
            elif constrain == 'sliding':
                constraintlist = list()
            elif ',' in constrain:
                log.info(logid+'Calculating probs for constraint ' + constrain)
                constraintlist = constrain.split(',')
            else:
                log.error(logid+'Could not compute constraints from input')
                return

            for fa in SeqIO.parse(seq,'fasta'):
                goi, chrom, strand = idfromfa(fa.id)

                outdict[goi]=manager.list()  # add a managed list to append results to

                if pattern and pattern not in goi:
                    next
                else:
                    log.info(logid+'Working on ' + goi)
                    xvals = []
                    for y in range(1,len(fa.seq)+1):
                        xvals.append(y)
                    xs = np.array(xvals)
                    xvals = []

                try:
                   if goi in constraintlist:
                       conslist = constraintlist[goi]
                   elif 'generic' in constraintlist:
                       conslist = constraintlist['generic']
                   elif 'NOCONS' in constraintlist:
                       conslist = constraintlist
                   elif constrain == 'sliding':
                       for start in range(1,len(fa.seq)-conslength+2):
                           end = start+conslength-1
                           conslist.append(str(start)+'-'+str(end))
                   else:
                       log.info(logid+'Could not map constraints to sequence for '+goi)
                       continue
                except Exception as err:
                        exc_type, exc_value, exc_tb = sys.exc_info()
                        tbe = tb.TracebackException(
                            exc_type, exc_value, exc_tb,
                        )
                        log.error(logid+''.join(tbe.format()))

                for entry in conslist:
                    if entry == 'NOCONS': # in case we just want to fold the sequence without constraints at all
                        gibbs_uc = [pool.apply_async(fold_unconstraint, args=(fa.seq))]
                        return (gibbs_uc)

                    else:
                        # we now have a list of constraints and for the raw seq comparison we only need to fold windows around these constraints
                        fstart, fend = [None,None]
                        if constrain == 'paired':
                            log.info(logid+'Calculating constraint pair\t' + entry)
                            fstart, fend, start, end = map(int,entry.split(',',1).split('-',1))
                        elif any(x in constrain for x in ['paired', 'Paired']):
                            fstart, fend, start, end = map(int,entry.split('-'))
                        else:
                            log.info(logid+'Calculating constraint\t' + entry)
                            start,end = map(int,entry.split('-',1))
                            #for all constraints we now extract subsequences to compare against
                            #we no longer fold the whole raw sequence but only the constraint region +- window size

                        if fstart and fend:
                            cons = str(start)+'-'+str(end)+':'+str(fstart)+'-'+str(fend)
                            const = np.sort(np.array([start, end, fstart, fend]))  # Is sorting ok?
                        else:
                            cons = str(start)+'-'+str(end)
                            const = np.array([start, end])

                        pool.apply_async(constrain_seq, args=(fa, start, end, conslength, const, cons, window, span, xs, unconstraint, paired, unpaired, save, outdir, genecoords, plot, outdict[goi]))

            pool.close()
            pool.join()

            for result in [outdict[gene] for gene in outdict]:
                for r in result:
                    write_out(r)

        log.info(logid+"DONE: output in: " + str(outdir))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

##### Functions #####

def constrain_seq(fa, start, end, conslength, const, cons, window, span, xs, unconstraint, paired, unpaired, save, outdir, genecoords, plot, outlist):
    #   DEBUGGING
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)

    logid = scriptname+'.constrain_seq: '
    try:
        goi, chrom, strand = idfromfa(fa.id)
        dist, fstart, fend = [None, None, None]
        if ':' in cons:
            start, end, fstart, fend = const
            dist=abs(start-fend)
            if dist > window:
                return(log.warning(logid+'Window '+str(window)+' too small for constraint distance '+str(dist)))
        else:
            start, end = const

        if dist:
            tostart = start - (window - abs(dist))
            toend = fend + (window - abs(dist)) + 1
            log.debug(logid+'Considering distance '+str(dist)+' between constraints '+str(const)+' for window extraction from '+str(tostart)+' to '+str(toend))
        else:
            tostart = start - window
            toend = end + window + 1

        #we no longer fold the whole sequence but only the constraint region +- window size
        if window is not None:
            tostart = start - window
            toend = end + window + 1
            if tostart < 0:
                tostart = 0
            if toend > len(fa.seq):
                toend = len(fa.seq)
            if fend is not None and toend < fend:
                log.error(logid+'Constraint '+str(cons)+' out of sequence range '+str(toend)+'!Skipping!')  # One of the constraints is outside the sequence window
                return
            seqtofold = str(fa.seq[tostart:toend+1])
        else:
            tostart = 0
            toend = 0
            seqtofold = str(fa.seq[:])

        if len(seqtofold) < window:
            log.error(logid+'Sequence of '+goi+' to short, seqlenght '+str(len(seqtofold))+' with window size '+str(window)+'!Skipping! ')
            return(logid+'Sequence of '+goi+' to short, seqlenght '+str(len(seqtofold))+' with window size '+str(window)+'!Skipping!')

        #check constraints
        if tostart < 0:
            tostart = 0
        if toend > len(fa.seq):
            toend = len(fa.seq)

        cstart = start-tostart
        cend =  end-tostart

        if cstart < 0 or cend > len(fa.seq):
            log.warning(logid+'start of constraint '+str(cstart)+' end of constraint '+str(cend)+' while length of sequence '+str(len(fa.seq))+'! Skipping!')
            return

        checklist = []

        checklist.append((cstart,cend))
        if fstart and fend:
            cfstart = fstart-tostart
            cfend = fend-tostart
            checklist.append((cfstart,cfend))
            checklist.append((cstart,cend,cfstart,cfend))

        #data
        data = {'seq' : seqtofold, 'stru' : []}

        #set model details
        md = RNA.md()
        md.max_bp_span = span

        log.debug(logid+'Constraints for '+goi+' are '+str(checklist))

        for check in checklist:
            s,e,os,oe = [None,None,None,None]
            gs, ge = map(int, genecoords[goi][0].split(sep='-'))
            if len(check) < 3:
                s,e = check
                sp, ep = [s+gs, e+gs]
                gtostart, gtoend = [tostart+gs, toend+gs]
                printcons = str.join('|',[str.join('-',[str(tostart), str(toend)]), str.join('-',[str(gtostart), str(gtoend)]), str.join('-',[str(s), str(e)]), str.join('-',[str(sp), str(ep)])])
            else:
                s,e,os,oe = check
                sp,ep,osp,oep = [s+gs, e+gs, os+gs, oe+gs]
                gtostart, gtoend = [tostart+gs, toend+gs]
                log.debug(logid+'PAIRED:'+';'.join(map(str,[s,e,os,oe,gs,ge,sp,ep,osp,oep,tostart,toend,gtostart,gtoend])))
                printcons = str.join('|',[str.join('-', [str(tostart), str(toend)]), str.join('-',[str(gtostart), str(gtoend)]), str.join(':',[str.join('-',[str(s),str(e)]), str.join('-',[str(os),str(oe)])]), str.join(':',[str.join('-',[str(sp),str(ep)]), str.join('-',[str(osp),str(oep)])])])

            # create new fold_compound object
            fc = RNA.fold_compound(data['seq'], md)
            fc_u = RNA.fold_compound(data['seq'], md)
            fc_p = RNA.fold_compound(data['seq'], md)

            # call pf and prop calculation
            gibbs = fc.pf()[1]
            bppm = get_bppm(fc.bpp(), cstart, cend)

            if bppm is None:
                log.error(logid+'Empty bpp matrix returned, stopping here!')
                #                sys.exit(logid+'Empty bpp matrix returned, stopping here!')
                return

            #enforce paired
            fc_p = constrain_paired(fc_p, s, e)
            #enforce unpaired
            fc_u = constrain_unpaired(fc_u, s, e)
            #calculate probs and nrg
            gibbs_u = calc_gibbs(fc_u)
            bppm_u = get_bppm(fc_u.bpp(), cstart, cend)
            dg_u = gibbs_u - gibbs

            gibbs_p = calc_gibbs(fc_p)
            bppm_p = get_bppm(fc_p.bpp(), cstart, cend)
            dg_p = gibbs_p - gibbs

            ###  access_energy([a, b]) = -RT log(prob([a, b]))
            ###  prob([a, b]) = sum_{s in S[a, b]} exp(-E(s)/RT) / sum_{s in S0} exp(-E(s)/RT)
            bpp = calc_bpp(bppm)
            bpp_u = calc_bpp(bppm_u)
            bpp_p = calc_bpp(bppm_p)
            nrg = calc_nrg(bpp)
            nrg_u = calc_nrg(bpp_u)
            nrg_p = calc_nrg(bpp_p)

            fn = 'constraint'

            if os and oe:
                fn = 'pairedconstraint'

            outlist.append([fa, fn, gibbs, '0', nrg, printcons, str(window), str(span), outdir, 'unconstraint'])  # unconstraint
            outlist.append([fa, fn, gibbs_u, dg_u, nrg_u, printcons, str(window), str(span), outdir, 'constraint_unpaired'])  # constraint_paired
            outlist.append([fa, fn, gibbs_p, dg_p, nrg_p, printcons, str(window), str(span), outdir, 'constraint_paired'])  # constraint_unpaired

            if os and oe:
                log.debug(logid+'Second constraint: '+str(','.join(map(str,[goi,data['seq'],len(data['seq']),s,e,os,oe]))))
                #enforce both constraints
                #enforce both paired
                fc_p = constrain_paired(fc_p, os, oe)
                #enforce both unpaired
                fc_u = constrain_unpaired(fc_u, os, oe)
                #calculate probs and nrg
                gibbs_u = calc_gibbs(fc_u)
                bppm_u = get_bppm(fc_u.bpp(), cstart, cend)
                dg_u = gibbs_u - gibbs

                gibbs_p = calc_gibbs(fc_p)
                bppm_p = get_bppm(fc_p.bpp(), cstart, cend)
                dg_p = gibbs_p - gibbs

                bpp = calc_bpp(bppm)
                bpp_u = calc_bpp(bppm_u)
                bpp_p = calc_bpp(bppm_p)
                nrg = calc_nrg(bpp)
                nrg_u = calc_nrg(bpp_u)
                nrg_p = calc_nrg(bpp_p)

                outlist.append([fa, fn, gibbs_u, dg_u, nrg_u, printcons, str(window), str(span), outdir, 'bothconstraint_unpaired'])  # bothconstraint_unpaired
                outlist.append([fa, fn, gibbs_p, dg_p, nrg_p, printcons, str(window), str(span), outdir, 'bothconstraint_paired'])  # bothconstraint_paired

                #enforce second constraint
                # First clear old constraints
                fc_u = RNA.fold_compound(data['seq'], md)
                fc_p = RNA.fold_compound(data['seq'], md)
                #enforce second paired
                fc_p = constrain_paired(fc_p, os, oe)
                #enforce second unpaired
                fc_u = constrain_unpaired(fc_u, os, oe)
                #calculate probs and nrg
                gibbs_u = calc_gibbs(fc_u)
                bppm_u = get_bppm(fc_u.bpp(), cstart, cend)
                dg_u = gibbs_u - gibbs

                gibbs_p = calc_gibbs(fc_p)
                bppm_p = get_bppm(fc_p.bpp(), cstart, cend)
                dg_p = gibbs_p - gibbs

                bpp = calc_bpp(bppm)
                bpp_u = calc_bpp(bppm_u)
                bpp_p = calc_bpp(bppm_p)
                nrg = calc_nrg(bpp)
                nrg_u = calc_nrg(bpp_u)
                nrg_p = calc_nrg(bpp_p)

                outlist.append([fa, fn, gibbs_u, dg_u, nrg_u, printcons, str(window), str(span), outdir, 'secondconstraint_unpaired'])  # secondconstraint_unpaired
                outlist.append([fa, fn, gibbs_p, dg_p, nrg_p, printcons, str(window), str(span), outdir, 'secondconstraint_paired'])  # secondconstraint_paired


        return outlist

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def constrain_temp(fa, temp, window, span, an, animations, xs, save, outdir, plot):
#   print('FOLDING ' + str(fa.seq) + ' at temp ' + str(temp))
#refresh model details
    logid = scriptname+'.constrain_temp: '
    log.info(logid+'Constraining Temp to ' + temp)
    try:
	    md = RNA.md()
	    md.max_bp_span = span
	    md.window_size = window
	    #set temperature
	    md.temperature = temp
	    #create new fold_compound objects
	    fc_t = RNA.fold_compound(str(fa.seq), md, RNA.OPTION_WINDOW)
	    #data
	    data = {'seq' : str(fa.seq), 'stru' : '', 'nrg' : ''}
	    #set model details
	    md = RNA.md()
	    # create new fold_compound object
	    fc = RNA.fold_compound(data['seq'], md, RNA.OPTION_PF)
	    gibbs = fc.pf()
	    bppm = fc.bpp()

	    data['stru'] = gibbs[0]
	    data['nrg'] = gibbs[1]

	    # call prop window calculation
	    log.debug(logid+';'.join(map(str,[gibbs, constrain, conslength])))

	    for item in bppm:
	        for i in range(int(constrain),int(constrain)+conslength):
	            log.debug(bppm.index(item), i, item[i])

	    at = up_to_array(data_t['up'], None, len(fa.seq))

	    diff_nt = an - at

	    #####Set temp, and rewrite save and plot for temp
	    if save:
	        log.info(logid+'SAVINGTEMP')
	        write_temp(fa, str(temp), data_t, diff_nt, str(window), outdir)

	    if plot == 'svg' or plot == 'png':
	        log.info(logid+'PLOTTING '+str(plot))
	        plot_temp(fa, at, temp, xs, plot, outdir)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def foldaround(seq, fc, pos, clength, gibbs, nrg):
#here we take the already constraint fc and constrain regions of length clength around it to see what happens at the original binding site

    logid = scriptname+'.foldaround: '
    try:
        cstart = pos
        cend = pos+clength-1
        fc = constrain_unpaired(fc, cstart, cend)
        gibbs_u = calc_gibbs(fc)
        bppm = get_bppm(fc.bpp(), cstart, cend)
        bpp = calc_bpp(bppm)
        ddg = gibbs_u - gibbs
        nrg_u = calc_nrg(bpp)
        nrg_diff = nrg_u - nrg

        return [gibbs_u, ddg, nrg_diff]

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def fold_unconstraint(seq):
    logid = scriptname+'.fold_unconstraint: '
    try:
        #set model details
        md = RNA.md()
        # create new fold_compound object
        fc = RNA.fold_compound(seq, md, RNA.OPTION_PF)
        # call prop window calculation
        gibbs_uc = fc.pf()[1]

        return(gibbs_uc)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def write_out(result):

    fa, fname, gibbs, ddg, nrg, const, window, span, outdir, condition = result
    logid = scriptname+'.write_out: '
    goi, chrom, strand = idfromfa(fa.id)
    temp_outdir = os.path.join(outdir,goi)
    try:
        if fname != 'STDOUT':
            if not os.path.exists(temp_outdir):
                os.makedirs(temp_outdir)
            if not os.path.exists(os.path.join(temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span])+'.gz')):
                o = gzip.open(os.path.join(temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span])+'.gz'), 'wb')
                o.write(bytes(str.join('\t',['Condition','FreeNRG(gibbs)','deltaG','OpeningNRG'])+'\n',encoding='UTF-8'))
                o.write(bytes(str.join('\t',[condition, str(gibbs), str(ddg), str(nrg), str(const)])+'\n',encoding='UTF-8'))
            else:
                log.info(os.path.join(temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span])+'.gz')+' exists, will append!')
                o = gzip.open(os.path.join(temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span])+'.gz'), 'ab')
                o.write(bytes(str.join('\t',[condition, str(gibbs), str(ddg), str(nrg), str(const)])+'\n',encoding='UTF-8'))
        else:
            print(str.join('\t',['FreeNRG(gibbs)','deltaG','OpeningNRG','Constraint']))
            print(str.join('\t',[str(gibbs), str(ddg), str(nrg), str(const)]))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        if args.loglevel != 'WARNING':
          streamlog = setup_multiprocess_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
          log = setup_multiprocess_logger(name=scriptname, log_file='logs/'+scriptname, logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)

        log.info(logid+'Running ConstraintFold on '+str(args.procs)+' cores')
        preprocess(args.sequence, args.window, args.span, args.unconstraint, args.unpaired, args.paired, args.length, args.gc, args.number, args.constrain, args.conslength, args.alphabet, args.plot, args.save, args.procs, args.vrna, args.temprange, args.outdir, args.genes, args.verbosity, args.pattern, args.cutoff)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
# ConstraintFold.py ends here
