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
import collections
import six
import functools

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
        exc_type,
        exc_value,
        exc_tb,
    )
    print("".join(tbe.format()), file=sys.stderr)


### Wrapper for Try/Except
def check_run(func):
    @functools.wraps(func)
    def func_wrapper(*args, **kwargs):
        logid = scriptn + ".Collection_func_wrapper: "
        try:
            return func(*args, **kwargs)
        except:                        
            return Exception
    return func_wrapper


### error_callback for pool.apply_async
def on_error(error):
    log.error(f"SOMETHING WENT WRONG, ERROR {type(error).__name__}('{error}') SHUTTING DOWN")
    raise error
    # queue.join()


##############################
####### VALIDITY CHECK #######
##############################


def isvalid(x=None):
    logid = scriptn + ".isvalid: "
    try:
        if x or x == 0:
            if x in ("None", "nan", "none", "NA", "NAN") or x is None or x is np.nan:
                return False
            else:
                return True
        else:
            return False
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def isinvalid(x=None):
    logid = scriptn + ".isinvalid: "
    try:
        if x or x == 0:
            if x in ("None", "nan", "none", "NA", "NAN") or x is None or x is np.nan:
                return True
            else:
                return False
        else:
            return True
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


##############################
####### Helper #######
##############################


def makeoutdir(outdir):
    logid = scriptn + ".makeoutdir: "
    try:
        if not os.path.isabs(outdir):
            outdir = os.path.abspath(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        return outdir
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def get_location(entry):
    logid = scriptn + ".get_location: "
    try:
        ret = list()
        start = end = strand = None
        start, end = map(int, entry.split(sep="|")[0].split(sep="-"))
        strand = str(entry.split(sep="|")[1])
        ret.extend([start, end, strand])

        if any([x == None for x in ret]):
            log.warning(logid + "Undefined variable: " + str(ret))

        log.debug(logid + str.join(" ", [str(entry), str(ret)]))
        return ret

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def expand_window(start, end, window, multiplyer, seqlen):
    logid = scriptn + ".expand_window: "
    try:
        tostart = start - multiplyer * window
        if tostart < 1:
            tostart = 1
        toend = end + multiplyer * window
        if toend > seqlen:
            toend = seqlen
        return [tostart, toend]
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def localize_window(start, end, window, seqlen):
    logid = scriptn + ".localize_window: "
    try:
        diff = start - window
        if diff < 1:
            locws = 1
        else:
            locws = diff

        locwe = (
            diff + 2 * window + (end - start)
        )  # this makes sure that if the start was trimmed, we do not just extend too much

        if locwe > seqlen:
            locwe = seqlen
        return [locws, locwe]
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


### Code utils


def print_globaldicts():
    logid = scriptn + ".print_globaldicts: "
    try:
        for name, value in globals().copy().items():
            print(name, value)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def print_globallists():
    logid = scriptn + ".print_globallists: "
    try:
        for name, value in globals().deepcopy().items():
            print(name, value)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


## DataStructure Helper
def merge_dicts(d, u):
    # https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    # python 3.8+ compatibility
    try:
        collectionsAbc = collections.abc
    except:
        collectionsAbc = collections

    for k, v in six.iteritems(u):
        dv = d.get(k, {})
        if not isinstance(dv, collectionsAbc.Mapping):
            d[k] = v
        elif isinstance(v, collectionsAbc.Mapping):
            d[k] = merge_dicts(dv, v)
        else:
            d[k] = v
    return d


# Code utils

# def print_globaldicts():
#     logid = scriptn+'.print_globaldicts: '
#     try:
#         for name, value in globals().copy().items():
#             print(name, value)
#     except Exception:
#         exc_type, exc_value, exc_tb = sys.exc_info()
#         tbe = tb.TracebackException(
#             exc_type, exc_value, exc_tb,
#             )
#         log.error(logid+''.join(tbe.format()))
#
#
# def print_globallists():
#     logid = scriptn+'.print_globallists: '
#     try:
#         for name, value in globals().deepcopy().items():
#             print(name, value)
#     except Exception:
#         exc_type, exc_value, exc_tb = sys.exc_info()
#         tbe = tb.TracebackException(
#             exc_type, exc_value, exc_tb,
#             )
#         log.error(logid+''.join(tbe.format()))

#
# Collection.py ends here
