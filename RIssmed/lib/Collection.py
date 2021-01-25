### Collection.py ---
##
## Filename: Collection.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Sep  6 09:02:18 2018 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Wed Dec 16 13:19:45 2020 (+0100)
##           By: Joerg Fallmann
##     Update #: 384
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
import logging
## other modules
import traceback as tb
import argparse
import numpy as np

############################################################
######################## FUNCTIONS #########################
############################################################

try:
    log = logging.getLogger(__name__)  # use module name
    scriptn = __name__  # os.path.basename(inspect.stack()[-1].filename).replace('.py', '')
except Exception:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    print(''.join(tbe.format()),file=sys.stderr)


##############################
########## ARGPARSE ##########
##############################

def parseargs_plcons():
    parser = argparse.ArgumentParser(description='Calculate base pairing probs of given seqs or random seqs for given window size, span and region.')
    parser.add_argument("-s", "--sequence", type=str, help='Sequence to fold')
    parser.add_argument("-w", "--window", type=int, default=240, help='Size of window')
    parser.add_argument("-l", "--span", type=int, default=60, help='Length of bp span')
    parser.add_argument("-u", "--region", type=int, default=1, help='Length of region')
    parser.add_argument("-m", "--multi", type=int, default=2, help='Multiplyer for window expansion')
    parser.add_argument("-c", "--cutoff", type=float, default=-0.01, help='Only print prob greater cutoff')
    parser.add_argument("-r", "--unconstraint", type=str, default='STDOUT', help='Print output of unconstraint folding to file with this name')
    parser.add_argument("-n", "--unpaired", type=str, default='STDOUT', help='Print output of unpaired folding to file with this name')
    parser.add_argument("-p", "--paired", type=str, default='STDOUT', help='Print output of paired folding to file with this name')
    parser.add_argument("-e", "--length", type=int, default=100, help='Length of randseq')
    parser.add_argument("--gc", type=int, default=0, help='GC content, needs to be %%2==0 or will be rounded')
    parser.add_argument("-b", "--number", type=int, default=1, help='Number of random seqs to generate')
    parser.add_argument("-x", "--constrain", type=str, default='sliding', help='Region to constrain, either sliding window (default) or region to constrain (e.g. 1-10) or path to file containing regions following the naming pattern $fastaID_constraints e.g. Sequence1_constraints, if paired, the first entry of the file will become a fixed constraint and paired with all the others, choices = [off,sliding,temperature, tempprobe, the string file or a filename, paired, or simply 1-10,2-11 or 1-10;15-20,2-11;16-21 for paired or ono(oneonone),filename to use one line of constraint file for one sequence from fasta]')
    parser.add_argument("-y", "--conslength", type=int, default=1, help='Length of region to constrain for slidingwindow')
    parser.add_argument("-t", "--temprange", type=str, default='', help='Temperature range for structure prediction (e.g. 37-60)')
    parser.add_argument("-a", "--alphabet", type=str, default='AUCG', help='alphabet for random seqs')
    #parser.add_argument("--plot", type=str, default='0', choices=['0','svg', 'png'], help='Create image of the (un-)constraint sequence, you can select the file format here (svg,png). These images can later on be animated with ImageMagick like `convert -delay 120 -loop 0 *.svg animated.gif`.')
    parser.add_argument("--save", type=int, default=0, choices=[0, 1], help='Save the output as numpy files only [0] or also as gzipped text files [1]')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory to write to')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processes to run this job with')
    parser.add_argument("--vrna", type=str, default='', help="Append path to vrna RNA module to sys.path")
    parser.add_argument("--pattern", type=str, default='', help="Helper var, only used if called from other prog where a pattern for files is defined")
    parser.add_argument("-g", "--genes", type=str, default='', help='Genomic coordinates bed for genes, either standard bed format or AnnotateBed.pl format')
    parser.add_argument("-v", "--verbosity", type=int, default=0, choices=[0, 1], help="increase output verbosity")
    parser.add_argument("--loglevel", type=str, default='WARNING', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")
    parser.add_argument("--logdir", type=str, default='LOGS', help="Set log directory")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def parseargs_collectpl():
    parser = argparse.ArgumentParser(description='Calculate the regions with highest accessibility diff for given Sequence Pattern')
    parser.add_argument("-p", "--pattern", type=str, default='250,150', help='Pattern for files and window, e.g. Seq1_30,250')
    parser.add_argument("-c", "--cutoff", type=float, default=.2, help='Cutoff for the definition of pairedness, if set to e.g. 0.2 it will mark all regions with probability of being unpaired >= cutoff as unpaired')
    parser.add_argument("-b", "--border", type=float, default=0.0, help='Cutoff for the minimum change between unconstraint and constraint structure, regions below this cutoff will not be returned as list of regions with most impact on structure.')
    parser.add_argument("-u", "--ulimit", type=int, default=1, help='Stretch of nucleotides used during plfold run (-u option)')
    parser.add_argument("-r", "--roi", type=str, default=None, help='Define Region of Interest that will be compared')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory to write to')
    parser.add_argument("-d", "--dir", type=str, default='', help='Directory to read from')
    parser.add_argument("-g", "--genes", type=str, help='Genomic coordinates bed for genes, either standard bed format or AnnotateBed.pl format')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processes to run this job with, only important of no border is given and we need to fold')
    parser.add_argument("--loglevel", type=str, default='WARNING', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")
    parser.add_argument("--logdir", type=str, default='LOGS', help="Set log directory")
    parser.add_argument("-w", "--padding", type=int, default=1, help='Padding around constraint that will be excluded from report, default is 1, so directly overlapping effects will be ignored')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def parseargs_foldcons():
    parser = argparse.ArgumentParser(description='Calculate base pairing probs of given seqs or random seqs for given window size, span and region.')
    parser.add_argument("-s", "--sequence", type=str, help='Sequence to fold')
    parser.add_argument("-w", "--window", type=int, default=None, help='Size of window')
    parser.add_argument("-l", "--span", type=int, default=None, help='Maximum bp span')
    parser.add_argument("-c", "--cutoff", type=float, default=-0.01, help='Only print prob greater cutoff')
    parser.add_argument("-r", "--unconstraint", type=str, default='STDOUT', help='Print output of unconstraint folding to file with this name, use "STDOUT" to print to STDOUT (default)')# or "add" to append to paired/unpaired output')
    parser.add_argument("-n", "--unpaired", type=str, default='STDOUT', help='Print output of unpaired folding to file with this name')
    parser.add_argument("-p", "--paired", type=str, default='STDOUT', help='Print output of paired folding to file with this name')
    parser.add_argument("-e", "--length", type=int, default=100, help='Length of randseq')
    parser.add_argument("--gc", type=int, default=0, help='GC content, needs to be %%2==0 or will be rounded')
    parser.add_argument("-b", "--number", type=int, default=1, help='Number of random seqs to generate')
    parser.add_argument("-x", "--constrain", type=str, default='sliding', help='Region to constrain, either sliding window (default) or region to constrain (e.g. 1-10) or path to file containing regions following the naming pattern $fastaID_constraints, if paired, the first entry of the file will become a fixed constraint and paired with all the others, e.g. Sequence1_constraints, choices = [off,sliding,temperature, tempprobe, file, paired, or simply 1-10,2-11 or 1-10;15-20,2-11:16-21 for paired]')
    parser.add_argument("-y", "--conslength", type=int, default=0, help='Length of region to constrain for slidingwindow')
    parser.add_argument("-t", "--temprange", type=str, default='', help='Temperature range for structure prediction (e.g. 37-60)')
    parser.add_argument("-a", "--alphabet", type=str, default='AUCG', help='alphabet for random seqs')
    #parser.add_argument("--plot", type=str, default='0', choices=['0','svg', 'png'], help='Create image of the (un-)constraint sequence, you can select the file format here (svg,png). These images can later on be animated with ImageMagick like `convert -delay 120 -loop 0 *.svg animated.gif`.')
    parser.add_argument("--save", type=int, default=1, help='Save the output as gz files')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory to write to')
    parser.add_argument("-z", "--procs", type=int, default=1, help='Number of parallel processed to run this job with')
    parser.add_argument("--vrna", type=str, default='', help="Append path to vrna RNA module to sys.path")
    parser.add_argument("--pattern", type=str, default='', help="Helper var, only used if called from other prog where a pattern for files is defined")
    parser.add_argument("-g", "--genes", type=str, default='', help='Genomic coordinates bed for genes, either standard bed format or AnnotateBed.pl format')
    parser.add_argument("-v", "--verbosity", type=int, default=0, choices=[0, 1], help="increase output verbosity")
    parser.add_argument("--loglevel", type=str, default='WARNING', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def parseargs_collectWindow():
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


##############################
####### VALIDITY CHECK #######
##############################

def isvalid(x=None):
    logid = scriptn+'.isvalid: '
    try:
        if x:
            if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
                return False
            else:
                return True
        else:
            return False
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def isinvalid(x=None):
    logid = scriptn+'.isinvalid: '
    try:
        if x:
            if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
                return True
            else:
                return False
        else:
            return True
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

##############################
####### Helper #######
##############################

def makeoutdir(outdir):
    logid = scriptn+'.makeoutdir: '
    try:
        if not os.path.isabs(outdir):
            outdir =  os.path.abspath(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        return outdir
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def get_location(entry):
    logid = scriptn+'.get_location: '
    try:
        ret = list()
        start = end = strand = None
        start, end = map(int, entry.split(sep='|')[0].split(sep='-'))
        strand = str(entry.split(sep='|')[1])
        ret.extend([start, end, strand])

        if any([x == None for x in ret]):
            log.warning(logid+'Undefined variable: '+str(ret))

        log.debug(logid+str.join(' ',[str(entry),str(ret)]))
        return ret

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def expand_window(start, end, window, multiplyer, seqlen):
    logid = scriptn+'.expand_window: '
    try:
        tostart = start - multiplyer*window
        if tostart < 1:
            tostart = 1
        toend = end + multiplyer*window
        if toend > seqlen:
            toend = seqlen
        return [tostart, toend]
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


def localize_window(start, end, window, seqlen):
    logid = scriptn+'.localize_window: '
    try:
        diff = start - window
        if diff < 1:
            locws = 1
        else:
            locws = diff

        locwe = diff + 2 * window + (end - start) # this makes sure that if the start was trimmed, we do not just extend too much

        if locwe > seqlen:
            locwe = seqlen
        return [locws, locwe]
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


### Code utils

def print_globaldicts():
    logid = scriptn+'.print_globaldicts: '
    try:
        for name, value in globals().copy().items():
            print(name, value)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


def print_globallists():
    logid = scriptn+'.print_globallists: '
    try:
        for name, value in globals().deepcopy().items():
            print(name, value)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

#
# Collection.py ends here
