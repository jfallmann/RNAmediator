#!/usr/bin/env python3
# Randseq.py ---
#
# Filename: Randseq.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Jul 11 13:29:38 2017 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Sep  1 10:14:44 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 174
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

from random import choice, shuffle
from io import StringIO
import gzip
# own
from RIssmed.RNAtweaks.RIssmedArgparsers import *
from RIssmed.RNAtweaks.logger import makelogdir, setup_logger
# Create log dir
makelogdir('LOGS')
# Define loggers
scriptname=os.path.basename(__file__)

def parseargs():
    parser = argparse.ArgumentParser(description='Generate random sequences of length l, if needed with gc content of g.')
    parser.add_argument("-l", "--length", type=int, default=100, help='Length of randseq')
    parser.add_argument("-g", "--gc", type=int, help='GC content, needs to be %2==0 or will be rounded')
    parser.add_argument("-n", "--number", type=int, default=1, help='Number of random seqs to generate')
    parser.add_argument("-o", "--outfile", type=str, default='Random', help='Name of output file for random sequences')
    parser.add_argument("-a", "--alphabet", type=str, default='AUCG', help='alphabet for random seqs')
    parser.add_argument("-v", "--verbosity", type=int, default=0, choices=[0, 1, 2], help="increase output verbosity")
    return parser.parse_args()

def createrandseq(length, gc, number, alphabet, outfile='Random', verbosity=False):
    try:
        nucs = list(alphabet)
        seqs=[]
        for i in range(number):
            if gc:
                content = (gc+(gc%2))/2
                rest    = (100-2*content)/2
                occ     = [];
                for x in nucs:
                    if x == 'C' or x == 'G':
                        occ.append(content*length/100)
                    else:
                        occ.append(rest*length/100)
                probs = iter(occ)
                seq = weightedrandseq(nucs, probs, length)
                header = ">Seq{i}_{gc}:random:nochrom:(.)\n".format(i=i+1,gc=gc)

            else:
                seq = randseq(alphabet, length)
                header = ">Seq{i}:random:nochrom:(.)\n".format(i=i+1)

            final = (''.join(seq))
            seqs.append(str("{header}{final}".format(header=header, final=final)))
        return seqs
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def randseq(items, length):
    try:
        l=''
        for i in range(length):
            l += str(''.join(choice(items)))
        return str(l)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def weightedrandseq(items, probs , length):
    try:
        seq=[]
        for i in items:
            weight=int(next(probs))
            for a in range(weight):
                seq+=i

        if len(seq) < length:
            for i in range(length-len(seq)):
                seq += str(''.join(choice(items)))

        shuffle(seq)
        return seq
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))


#choices(items,weights=w,k=nr)

if __name__ == '__main__':
    try:
        args=parseargs()
        log = setup_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='WARNING')
        rand = "\n".join(createrandseq(args.length, args.gc, args.number, args.alphabet, args.outfile, args.verbosity))
        seq = StringIO(rand)
        o=gzip.open(args.outfile+'.fa.gz','wb')
        o.write(bytes(rand,encoding='UTF-8'))
        o.close()
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))
#
# Randseq.py ends here
