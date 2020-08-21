# Random.py ---
#
# Filename: Random.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Fri Aug 21 09:40:58 2020 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Fri Aug 21 10:13:59 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 1
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
from lib.logger import *
##other modules
import numpy as np
import heapq
from operator import itemgetter
from natsort import natsorted, ns
import traceback as tb
from io import StringIO
#Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq

################
# Random Seqs  #
################

def create_kmers(choices, length):
    logid = scriptn+'.create_kmers: '
    try:
        #     choices=['A','T','G','C']
        bases = list(choices)
        k = length

        return [''.join(p) for p in itertools.product(bases, repeat=k)]

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def randseq(alphabet, length):
    logid = scriptn+'.randseq: '
    try:
        l=''
        for i in range(length):
            l += str(''.join(choice(alphabet)))
            return str(l)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def weightedrandseq(alphabet, probs , length):
    logid = scriptn+'.weightedrandseq: '
    try:
        seq=[]
        for i in alphabet:
            weight=int(next(probs))
            for a in range(weight):
                seq+=i

        if len(seq) < length:
            for i in range(length-len(seq)):
                seq += str(''.join(choice(alphabet)))

        shuffle(seq)
        return seq
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


#
# Random.py ends here
