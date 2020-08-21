# NPtweaks.py ---
#
# Filename: NPtweaks.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Fri Aug 21 10:20:20 2020 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Fri Aug 21 10:21:13 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 2
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
import traceback as tb
import numpy as np

####################
# Numpy processing #
####################

def toarray(file, ulim=None):
    logid = scriptn+'.toarray: '
    try:
        if not ulim:
            ulim = 1
        x = np.loadtxt(str(file), usecols = (ulim), delimiter = '\t', unpack = True, converters = {ulim: lambda s: convertcol(s.decode("utf-8"))})
        return x
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def convertcol(entry):
    logid = scriptn+'.convertcol: '
    try:
        if isinvalid(entry):
#       if entry is None or entry == 'NA' or entry == 'nan' or entry is np.nan:
            return np.nan
        else:
            return float(entry)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


#
# NPtweaks.py ends here
