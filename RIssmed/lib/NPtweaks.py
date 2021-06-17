# NPtweaks.py ---
#
# Filename: NPtweaks.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Fri Aug 21 10:20:20 2020 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Sep  1 10:30:36 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 6
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
import os
import sys
import inspect
## other modules
import traceback as tb
import numpy as np
# own
import logging
from RIssmed.RNAtweaks.RIssmedArgparsers import *

####################
# Numpy processing #
####################

try:
    log = logging.getLogger(__name__)  # use module name
    scriptn = os.path.basename(inspect.stack()[-1].filename).replace('.py', '')
    log.debug('LOGGING IN NPtweaks'+str(scriptn)+str(log)+str(log.handlers))
except Exception:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    print(''.join(tbe.format()),file=sys.stderr)


def toarray(file, ulim=None):
    logid = scriptn+'.toarray: '
    try:
        if not ulim:
            ulim = 1
        x = np.loadtxt(str(file), usecols = (ulim), delimiter = '\t', unpack = True, converters = {ulim: lambda s: convertcol(s.decode("utf-8"))})
        return x
    except Exception:
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
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


#
# NPtweaks.py ends here
