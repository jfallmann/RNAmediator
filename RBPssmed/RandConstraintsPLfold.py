#!/usr/bin/env python

# RandConstraintsPLfold.py ---
#
# Filename: RandConstraintsPLfold.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Fri Jun 30 11:18:54 2017 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Mon Sep 24 13:08:14 2018 (+0200)
#           By: Joerg Fallmann
#     Update #: 15
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

import argparse
import os, sys, inspect
import time
import regex as re
#load own modules
cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) )),"../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
from Collection import *

#sys.path.append("/usr/local/ViennaRNA-2.1.7")
#import RNA

def parseargs():
	parser = argparse.ArgumentParser(description='Fold Fasta, add random constraints, fold again, plot changes.')
	parser.add_argument("-f", "--fasta", type=str, default=False, action='store_true', help='Read sequence from file')
	parser.add_argument("-p", "--package", type=str, default='vrna', help='Path to ViennaRNA (default: vrna)')
	parser.add_argument("-n", "--number", type=int, default=100, help='Number of random constraint iterations to run')
	parser.add_argument("-l", "--length", type=int, default=8, help='Length of random constraint to add')
	parser.add_argument("-u", "--unpaired", type=int, default=1, action='store_true', help='Set unpaired span for PLfold')
	parser.add_argument("-w", "--window", type=int, default=240, action='store_true', help='Set window length for PLfold')
	parser.add_argument("-s", "--span", type=int, default=160, action='store_true', help='Set span length for PLfold')
    parser.add_argument("--vrna", type=str, default='', help="Append path to vrna RNA module to sys.path")
	args = parser.parse_args()

def randconst(number, length, fasta, unpaired, window, span, vrna)

	print("# Options: numberofiterations={0:d}, constraint={1:d}, fasta={2:}, unpaired={3:d}, window={4:d}, span={5:d}".format(number, length, fasta, unpaired, window, span))

	f = parse_fasta(fasta)

	for line in f:
		if (re.match(r'^>/',line)):
			name = line[1:]
		else:
			seq = line
			open(name+'.nrgs',a)
			# compute minimum free energy (mfe) and corresponding structure
			(bps) = RNA.plfold(sequence)
			nrg =

def parse_fasta(fasta):
	f = open(fasta,'r')
	return f


if __name__ == '__main__':
    args=parseargs()
    randconst(args.number, args.length, args.fasta, args.unpaired, args.window, args.span, args.vrna)
#
# RandConstraintsPLfold.py ends here
