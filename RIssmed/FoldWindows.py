#!/usr/bin/env python3
## FoldWindows.py ---
##
## Filename: FoldWindows.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Sep  6 09:02:18 2018 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Thu Jul 16 08:48:06 2020 (+0200)
##           By: Joerg Fallmann
##     Update #: 163
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

#other modules
import argparse
import pprint
from io import StringIO
import time
import math
import gzip
import importlib
import multiprocessing
from multiprocessing import get_context
import traceback as tb
import shutil
import shlex
#Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq
#numpy and matplolib and pyplot
import numpy as np
from random import choices, choice, shuffle # need this if tempprobing was choosen
from Randseq import createrandseq

#Logging
import datetime
from lib.logger import makelogdir, setup_multiprocess_logger

scriptname = os.path.basename(__file__).replace('.py','')
makelogdir('LOGS')
if not os.path.isfile(os.path.abspath('LOGS/'+scriptname+'.log')):
    open('LOGS/'+scriptname+'.log','a').close()
else:
    ts = str(datetime.datetime.fromtimestamp(os.path.getmtime(os.path.abspath('LOGS/'+scriptname+'.log'))).strftime("%Y%m%d_%H_%M_%S"))
    shutil.copy2('LOGS/'+scriptname+'.log','LOGS/'+scriptname+'_'+ts+'.log')

logfile = 'LOGS/'+scriptname+'.log'
log = setup_multiprocess_logger(name=scriptname, log_file=logfile, filemode='a', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M')
log = setup_multiprocess_logger(name='', log_file='stderr', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M')

##load own modules
from lib.Collection import *

def parseargs():
    parser = argparse.ArgumentParser(description='Calculate base pairing probs of given seqs or random seqs for given window size, span and region.')
    parser.add_argument("-s", "--sequence", type=str, help='Sequence to fold')
    parser.add_argument("-w", "--window", type=int, default=240, help='Size of window')
    parser.add_argument("-l", "--span", type=int, default=60, help='Length of bp span')
    parser.add_argument("-u", "--region", type=int, default=1, help='Length of region')
    parser.add_argument("-t", "--printto", type=str, default='STDOUT', help='Print output of folding to file with this name')
    parser.add_argument("-e", "--length", type=int, default=100, help='Length of randseq')
    parser.add_argument("-g", "--gc", type=int, default=0, help='GC content, needs to be %2==0 or will be rounded')
    parser.add_argument("-b", "--number", type=int, default=1, help='Number of random seqs to generate')
    parser.add_argument("-a", "--alphabet", type=str, default='AUCG', help='alphabet for random seqs')
    #parser.add_argument("--plot", type=str, default='0', choices=['0','svg', 'png'], help='Create image of the (un-)constraint sequence, you can select the file format here (svg,png). These images can later on be animated with ImageMagick like `convert -delay 120 -loop 0 *.svg animated.gif`.')
    parser.add_argument("--save", type=int, default=1, help='Save the output as gz files')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory to write to')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processed to run this job with')
    parser.add_argument("--vrna", type=str, default='', help="Append path to vrna RNA module to sys.path")
    parser.add_argument("--pattern", type=str, default='', help="Helper var, only used if called from other prog where a pattern for files is defined")
    parser.add_argument("-v", "--verbosity", type=int, default=0, choices=[0, 1], help="increase output verbosity")

    return parser.parse_args()

def fold(sequence, window, span, region, printto, length, gc, number, alphabet, save, procs, vrna, outdir, verbosity=False, pattern=None):

    logid = scriptname+'.fold: '
    #set path for output
    if outdir:
        printlog(outdir)
        if not os.path.isabs(outdir):
            outdir =  os.path.abspath(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    else:
        outdir = os.path.abspath(os.getcwd())
#set path for VRNA lib
    if vrna:
        sys.path=[vrna] + sys.path
    else:
        sys.path=["/scratch/fall/VRNA/249/lib/python3.6/site-packages"] + sys.path
    try:
        global RNA
        RNA = importlib.import_module('RNA')
        globals().update(
            {n: getattr(RNA, n) for n in RNA.__all__}
            if hasattr(RNA, '__all__')
            else {k: v for (k, v) in RNA.__dict__.items() if not k.startswith('_')
                  })
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

    if ( sequence == 'random' ):
        rand = "\n".join(createrandseq(length, gc, number, alphabet))
        seq = StringIO(rand)
        o = gzip.open('Random.fa.gz','wb')
        o.write(bytes(rand,encoding='UTF-8'))
        o.close()

    elif (os.path.isfile(sequence)):
        if '.gz' in sequence :
            seq = gzip.open(sequence,'rt')
        else:
            seq = open(sequence,'rt')
    else:
        header = ">Seq1:default:nochrom:(.)"
        s = sequence
        seq = StringIO("{header}\n{s}".format(header=header, s=s))

    for fa in SeqIO.parse(seq,'fasta'):
        try:
            goi, chrom = fa.id.split(':')[::2]
            strand = str(fa.id.split(':')[3].split('(')[1][0])
        except:
            printlog('Fasta header is not in expected format, you will loose information on strand and chromosome')
            try:
                goi = fa.id
                chrom, strand = ['na','na']
            except:
                printlog('Could not assign any value from fasta header, please check your fasta files')

        if pattern and pattern not in goi:
            next
        else:
            printlog('Working on ' + goi)
#set kT for nrg2prob and vice versa calcs
            kT = 0.61632077549999997
# Create process pool with processes
            num_processes = procs or 1
            with get_context("spawn").Pool(processes=num_processes-1, maxtasksperchild=1) as pool:

                try:
                    for reg in range(0,len(fa.seq)-window+1):
    #subseq
                        seqtofold = str(fa.seq[reg:reg+window])
                        pool.apply_async(fold_windows, args=(fa, seqtofold, reg, window, span, region, save, printto, outdir))
                except Exception:
                    exc_type, exc_value, exc_tb = sys.exc_info()
                    tbe = tb.TracebackException(
                        exc_type, exc_value, exc_tb,
                        )
                    log.error(logid+''.join(tbe.format()))

                pool.close()
                pool.join()

    printlog("DONE: output in: " + str(outdir))

##### Functions #####
def fold_windows(fa, seq, reg, window, span, region, save, printto, outdir):
#   DEBUGGING
#   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    logid = scriptname+'.fold_windows: '
    try:
        goi, chrom = fa.id.split(':')[::2]
        strand = str(fa.id.split(':')[3].split('(')[1][0])
    except:
        goi = fa.id
        chrom, strand = ['na','na']

    try:
#model details
        md = RNA.md()
        md.max_bp_span = span
        md.window_size = window

# create new fold_compound object
        fc = RNA.fold_compound(seq, md, RNA.OPTION_WINDOW)
# call prop window calculation
        data = { 'up': [] }
        fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data)
        if save:
            write_out(fa, printto, str(reg), len(seq), data, int(region), str(window), str(span), outdir)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def up_callback(v, v_size, i, maxsize, what, data):
    logid = scriptname+'.up_callback: '
    try:
        if what & RNA.PROBS_WINDOW_UP:
#        data['up'].extend([{ 'i': i, 'up': v}])
            data['up'].extend([v])
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def print_region_up(data=None, seqlength=None, region=None, winnr=None):
#   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
#   pp.pprint(data)
    logid = scriptname+'.print_region_up: '
    try:
        ups=''
        x = int(region)
        print('Printing region '+str(region)+' of data '+str(data[0]))
        for i in range(int(seqlength)):
            if data[i][x] is None or data[i][x] is np.nan or data[i][x] is 'nan':
                data[i][x] = 'NA'
            else:
                data[i][x] = round(data[i][x],7)
            ups+=str(i+1+winnr)+"\t"+str(data[i][x])+"\n"
        return ups
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def print_up(data=None, seqlength=None, region=None, winnr=None):
#   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    #   pp.pprint(data)
    logid = scriptname+'.print_up: '
    try:
        ups=''
        for i in range(int(seqlength)):
            for x in range(1,region+1):
                if data[i][x] is None or data[i][x] is np.nan or data[i][x] is 'nan':
                    data[i][x] = 'NA'
                else:
                    data[i][x] = round(data[i][x],7)
            ups+=str(i+1+int(winnr))+"\t"+"\t".join(map(str,data[i][1:region+1]))+"\n"
        return ups
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def up_to_array(data=None, region=None, seqlength=None):
#   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
#   pp.pprint(data[165553:165588])
    logid = scriptname+'.up_to_array: '
    try:
        entries=[]
        if not seqlength:
            seqlength = len(data)

        for i in range(seqlength):
            if data[i][region] is None or data[i][region] is 'NA' or data[i][region] is 'nan':
                data[i][region] = np.nan
                entries.append(data[i][region])
            else:
                entries.append(round(data[i][region],8))

        return np.array(entries)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def npprint(a, o=None):#, format_string ='{0:.2f}'):
    logid = scriptname+'.npprint: '
    try:
        out=''
        it = np.nditer(a, flags=['f_index'])
        while not it.finished:
            out += "%d\t%0.7f" % (it.index+1,it[0])+"\n"
            it.iternext()
        if o:
            o.write(bytes(out,encoding='UTF-8'))
        else:
            print(out)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def write_out(fa, printto, winnr, seqlen, data, region, window, span, outdir):
    logid = scriptname+'.write_out: '
    try:
        goi, chrom = fa.id.split(':')[::2]
        strand = str(fa.id.split(':')[3].split('(')[1][0])
    except:
        goi = fa.id
        chrom, strand = ['na','na']

    try:
        if printto != 'STDOUT':
            bakdir=os.path.abspath(os.getcwd())
            os.chdir(outdir)
            if not os.path.exists('Window_'+winnr+'_'+goi+'_'+chrom+'_'+strand+'_'+'_'+str(window)+'_'+str(span)+'.gz'):
                o = gzip.open(('Window_'+winnr+'_'+goi+'_'+chrom+'_'+strand+'_'+'_'+str(window)+'_'+str(span)+'.gz'), 'wb')
                o.write(bytes(print_up(data['up'],seqlen,region,winnr),encoding='UTF-8'))
            os.chdir(bakdir)
        else:
            print (print_up(data,seqlen,region))
        return 1
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def read_precalc_fold(data, name, fa):
    logid = scriptname+'.read_precalc_fold: '
    try:
        for i in range(len(fa.seq)):
            data.append([])
            data[i] = []
        with gzip.open(name,'rt') as o:
            for line in o:
                cells = line.rstrip().split('\t')
                data[int(cells[0])-1].append([])
                data[int(cells[0])-1][0] = None
                for a in range(1,len(cells)):
                    data[int(cells[0])-1].append([])
                    data[int(cells[0])-1][a] = float(cells[a])
        return data
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def checkexisting(fa, region, winnr, window, span, outdir):
    logid = scriptname+'.checkexisting: '
    try:
        goi, chrom = fa.id.split(':')[::2]
        strand = str(fa.id.split(':')[3].split('(')[1][0])
    except:
        goi = fa.id
        chrom, strand = ['na','na']
    try:
        if os.path.exists(os.path.join(outdir,'Window_'+winnr+'_'+goi+'_'+chrom+'_'+strand+'_'+'_'+str(window)+'_'+str(span)+'.gz')):
            return True
        else:
            return False
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

if __name__ == '__main__':

    try:
        scriptname=os.path.basename(__file__).replace('.py','')
        args=parseargs()
        logid = scriptname+'.main: '

        log.setLevel(args.loglevel)
        log.info(logid+'Running '+scriptname+' on '+str(args.procs)+' cores.')
        log.info(logid+'CLI: '+sys.argv[0]+'{}'.format(' '.join( [shlex.quote(s) for s in sys.argv[1:]] )))

        fold(args.sequence, args.window, args.span, args.region, args.printto, args.length, args.gc, args.number, args.alphabet, args.save, args.procs, args.vrna, args.outdir, args.verbosity, args.pattern)
            except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

# FoldWindows.py ends here
