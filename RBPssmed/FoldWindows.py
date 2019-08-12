#!/usr/bin/env python3
# ConstraintPLFold.py ---
#
# Filename: ConstraintPLFold.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Mon Oct 16 17:18:42 2017 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Mon Nov 12 19:16:00 2018 (+0100)
#           By: Joerg Fallmann
#     Update #: 2053
# URL:
# Doc URL:
# Keywords:
# Commentary:
# TODO:
# Animations per constraint if not sliding
# Constrain temp needs option to randomly select windows for long sequences
#

# TODO Log:
# Statistics about results
# Fold Windows
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

# IMPORTS
import argparse
import pprint
from io import StringIO
import time
import math
import os, sys, inspect
import gzip
import importlib
import multiprocessing
from multiprocessing import Manager
import traceback as tb
#Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq
#numpy and matplolib and pyplot
import numpy as np
import matplotlib
from matplotlib import animation #, rc
import matplotlib.pyplot as plt
from random import choices, choice, shuffle # need this if tempprobing was choosen
#load own modules
cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) )),"../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
from Collection import *
from Randseq import createrandseq

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
    parser.add_argument("--plot", type=str, default='0', choices=['0','svg', 'png'], help='Create image of the (un-)constraint sequence, you can select the file format here (svg,png). These images can later on be animated with ImageMagick like `convert -delay 120 -loop 0 *.svg animated.gif`.')
    parser.add_argument("--save", type=int, default=1, help='Save the output as gz files')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory to write to')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processed to run this job with')
    parser.add_argument("--vrna", type=str, default='', help="Append path to vrna RNA module to sys.path")
    parser.add_argument("--pattern", type=str, default='', help="Helper var, only used if called from other prog where a pattern for files is defined")
    parser.add_argument("-v", "--verbosity", type=int, default=0, choices=[0, 1], help="increase output verbosity")

    return parser.parse_args()

def fold(sequence, window, span, region, printto, length, gc, number, alphabet, plot, save, procs, vrna, outdir, verbosity=False, pattern=None):
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
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

    if ( plot == '0' and not save):
        raise ValueError('Neither plot nor save are active, this script will take a long time and produce nothing, please activate at least one of the two!')

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
##prepare plots
#set kT for nrg2prob and vice versa calcs
            kT = 0.61632077549999997
# Create process pool with processes
            num_processes = procs or 1
            pool = multiprocessing.Pool(processes=num_processes)

            try:
                for reg in range(0,len(fa.seq)-window+1):
#subseq
                    seqtofold = str(fa.seq[reg:reg+window])
                    pool.apply_async(fold_windows, args=(fa, seqtofold, reg, window, span, region, save, printto, outdir, plot))
            except Exception as err:
                exc_type, exc_value, exc_tb = sys.exc_info()
                tbe = tb.TracebackException(
                    exc_type, exc_value, exc_tb,
                    )
                with open('error','a') as h:
                    print(''.join(tbe.format()), file=h)

            pool.close()
            pool.join()

    printlog("DONE: output in: " + str(outdir))

##### Functions #####
def fold_windows(fa, seq, reg, window, span, region, save, printto, outdir, plot):
#   DEBUGGING
#   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)

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

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def up_callback(v, v_size, i, maxsize, what, data):
    if what & RNA.PROBS_WINDOW_UP:
#        data['up'].extend([{ 'i': i, 'up': v}])
        data['up'].extend([v])

def print_region_up(data=None, seqlength=None, region=None, winnr=None):
#   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
#   pp.pprint(data)
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
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def print_up(data=None, seqlength=None, region=None, winnr=None):
#   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    #   pp.pprint(data)
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
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def up_to_array(data=None, region=None, seqlength=None):
#   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
#   pp.pprint(data[165553:165588])
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
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def npprint(a, o=None):#, format_string ='{0:.2f}'):
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
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def write_out(fa, printto, winnr, seqlen, data, region, window, span, outdir):
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
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def read_precalc_fold(data, name, fa):
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
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def checkexisting(fa, region, winnr, window, span, outdir):
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
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

if __name__ == '__main__':
    args=parseargs()
    fold(args.sequence, args.window, args.span, args.region, args.printto, args.length, args.gc, args.number, args.alphabet, args.plot, args.save, args.procs, args.vrna, args.outdir, args.verbosity, args.pattern)
# ConstraintPLFold.py ends here
