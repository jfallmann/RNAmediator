# DStweaks.py ---
#
# Filename: DStweaks.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Fri Aug 21 10:18:26 2020 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Sep  1 10:30:50 2020 (+0200)
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

##other modules
import traceback as tb
import numpy as np
import heapq
from operator import itemgetter

# own
import logging

# from Tweaks.Collection import *

################
#  DS tweaker  #
################

try:
    log = logging.getLogger(__name__)  # use module name
    scriptn = os.path.basename(inspect.stack()[-1].filename).replace('.py', '')
    log.debug('LOGGING IN DStweaks' + str(scriptn) + str(log) + str(log.handlers))
except Exception:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type,
        exc_value,
        exc_tb,
    )
    print(''.join(tbe.format()), file=sys.stderr)


def removekey(d, key):
    logid = scriptn + '.removekey: '
    try:
        r = dict(d)
        del r[key]
        return r
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + ''.join(tbe.format()))


def getlowest_list(a, n):
    logid = scriptn + '.getlowest_list: '
    try:
        if n > len(a) - 1:
            b = len(a) - 1
        else:
            b = n
        if len(a) > 0 and n > 0:
            return list(np.partition(a, b)[:n])
        else:
            return list(None for i in range(n))
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + ''.join(tbe.format()))


def gethighest_list(a, n):
    logid = scriptn + '.gethighest_list: '
    try:
        if len(a) - n < 0:
            b = len(a) - 1
        else:
            b = len(a) - n
        if len(a) > 0 and n > 0:
            return list(np.partition(a, b)[-n:])
        else:
            return list(None for i in range(n))
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + ''.join(tbe.format()))


def getlowest_dict(a, n):
    logid = scriptn + '.getlowest_dict: '
    try:
        if n > len(a):
            b = len(a)
        else:
            b = n
        if len(a) > 0:
            return dict(heapq.nsmallest(b, a.items(), key=itemgetter(1)))
        else:
            return dict({i: None for i in range(n)})
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + ''.join(tbe.format()))


def gethighest_dict(a, n):
    logid = scriptn + '.gethighest_dict: '
    try:
        if n > len(a):
            b = len(a)
        else:
            b = n
        if len(a) > 0:
            return dict(heapq.nlargest(b, a.items(), key=itemgetter(1)))
        else:
            return dict({i: None for i in range(n)})
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + ''.join(tbe.format()))


#
# DStweaks.py ends here
