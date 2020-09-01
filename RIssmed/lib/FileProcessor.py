# FileProcessor.py ---
#
# Filename: FileProcessor.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Fri Aug 21 10:21:43 2020 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Sep  1 10:30:58 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 8
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

### IMPORTS
import os, sys, inspect
##other modules
import traceback as tb
import numpy as np
import heapq
from operator import itemgetter
from natsort import natsorted, ns
import re
import pprint
from io import StringIO
import gzip
import math
import datetime
from collections import defaultdict
#Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq
#own
from lib.logger import *
from lib.Collection import *

####################
# FILE processing  #
####################

try:
    scriptn = os.path.basename(inspect.stack()[-1].filename).replace('.py','')
    #scriptn = os.path.basename(__file__).replace('.py','')

except Exception as err:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    print(''.join(tbe.format()),file=sys.stderr)

def backup(file):
    logid = scriptn+'.backup: '
    try:
        if os.path.exists(file):
            os.rename(file,file+'.bak')
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

### SEQUENCE
def parseseq(sequence):
    logid = scriptn+'.parseseq: '
    try:
        if (isinstance(sequence, StringIO)):
            seq = sequence

        elif ( isinstance(sequence, str) and sequence == 'random' ):
            rand = "\n".join(createrandseq(length, gc, number, alphabet))
            seq = StringIO(rand)
            o = gzip.open('Random.fa.gz','wb')
            o.write(bytes(rand,encoding='UTF-8'))
            o.close()

        elif (isinstance(sequence, str) and os.path.isfile(sequence)):
            if '.gz' in sequence :
                seq = gzip.open(sequence,'rt')
            else:
                seq = open(sequence,'rt')
        else:
            header = ">Seq1:default:nochrom:(.)"
            s = sequence
            seq = StringIO("{header}\n{s}".format(header=header, s=s))

        return seq

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

### FASTA
def idfromfa(id):
    logid = scriptn+'.idfromfa: '
    goi, chrom, strand = [None, None, None]
    id = id.replace('_','-')
    try:
        goi, chrom = id.split(':')[::2]
        strand = str(id.split(':')[3].split('(')[1][0])
    except:
        log.error(logid+'Fasta header is not in expected format, you will loose information on strand and chromosome')
        goi = id
        chrom, strand = ['na','na']

    if goi and chrom and strand:
        return [str(goi), str(chrom), str(strand)]
    else:
        log.error(logid+'Could not assign any value from fasta header, please check your fasta files')
        sys.exit('Could not assign any value from fasta header, please check your fasta files')


###BED
def parse_annotation_bed(bed, annotated=None):
    logid = scriptn+'.parse_annotation_bed: '
    anno = defaultdict(list)
    if os.path.isfile(os.path.abspath(bed)):
        if '.gz' in bed:
            f = gzip.open(os.path.abspath(bed),'rt')
        else:
            f = open(os.path.abspath(bed),'rt')
    else:
        f = bed
    try:
        for line in f:
            entries = line.rstrip().split('\t')
            goi = entries[3]
            strand = entries[5]
            if annotated:
                start = int(entries[10])+1
                end   = int(entries[11])
                strand = entries[14]
            else:
                start = int(entries[1])+1
                end   = int(entries[2])
            anno[str(goi)].append('|'.join(['-'.join([str(start),str(end)]),strand]))  # Need strand info here!
        return anno
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def readConstraintsFromBed(bed, linewise=None):
    logid = scriptn+'.readConstraintsFromBed: '
    cons = defaultdict(list)
    try:
        for line in bed:
            entries = line.rstrip().split('\t')
            start = int(entries[1])+1
            end = entries[2]
            goi = entries[3]
            strand = entries[5]
            if linewise:
                cons['lw'].append('|'.join(['-'.join([str(start),str(end)]),strand]))
            else:
                cons[str(goi)].append('|'.join(['-'.join([str(start),str(end)]),strand]))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def readPairedConstraintsFromBed(bed, linewise=None):
    logid = scriptn+'.readPairedConstraintsFromBed: '
    cons = defaultdict(list)
    try:
        for line in bed:
            entries = line.rstrip().split('\t')
            if len(entries) % 2:
                raise Exception('Unbalanced paired bed, please make sure the paired bed consists of equal number of fields for both constraint entries')
            else:
                second = int((len(entries)/2)+1)
            if int(entries[1]) > -1 and int(entries[second]) > -1:
                start_one = int(entries[1])+1
                end_one = entries[2]
                goi = entries[3]
                strand = entries[5]
                start_two = int(entries[second])+1
                end_two = int(entries[second+1])
                if linewise:
                    cons['lw'].append(':'.join(['|'.join(['-'.join([str(start_one), str(end_one)]),strand]), '|'.join(['-'.join([str(start_two), str(end_two)]),strand])]))
                else:
                    cons[str(goi)].append(':'.join(['|'.join(['-'.join([str(start_one), str(end_one)]),strand]), '|'.join(['-'.join([str(start_two), str(end_two)]),strand])]))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

### CSV

def readConstraintsFromCSV(csv, linewise=None):
    logid = scriptn+'.readConstraintsCSV: '
    cons = defaultdict(
        lambda: defaultdict(list)
    )

    try:
        for line in csv:
            entries = split(',',line.rstrip())
            start = entries[1]
            end   = entries[2]
            strand = entries[5]
            if linewise:
                cons['def'].append('|'.join(['-'.join([str(start),str(end)]),strand]))
            else:
                cons[entries[3]].append('|'.join(['-'.join([str(start),str(end)]),strand]))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

### GENERIC

def readConstraintsFromGeneric(generic, linewise=None):
    logid = scriptn+'.readConstraintsFromGeneric: '
    cons = defaultdict(
        lambda: defaultdict(list)
    )

    try:
        for line in csv:
            entries = re.split(r'[ ,|;"]+', line.rstrip())
            if len(entries > 3):
                start = entries[2]
                end = entries[3]
                strand = entries[5]
                if linewise:
                    cons['lw'].append('|'.join(['-'.join([str(start),str(end)]),strand]))
                else:
                    cons[entries[0]].append('|'.join(['-'.join([str(start),str(end)]),strand]))
            else:
                start = entries[2]
                end = entries[3]
                strand = '.'
                if linewise:
                    cons['lw'].append('|'.join(['-'.join([str(start),str(end)]),strand]))
                else:
                    cons['generic'].append('|'.join(['-'.join([str(start),str(end)]),strand]))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))



#
# FileProcessor.py ends here
