## IMPORTS
import os
import sys
import multiprocessing
from multiprocessing import set_start_method, get_context, get_start_method
from io import StringIO
import gzip
import importlib
import traceback as tb
import shlex
# Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq
# numpy
import numpy as np
from random import choice  # need this if tempprobing was choosen
# RNA
import RNA
# Logging
import datetime
import logging
from lib.logger import makelogdir, makelogfile, listener_process, listener_configurer, worker_configurer
# load own modules
from lib.Collection import *
from lib.FileProcessor import *
from lib.RNAtweaks import *
from lib.NPtweaks import *
import errno


def run_api_cli_test():
    #dominik code here
    return None


def test_cli(seq, region, window, span, locws=None, locwe=None):
    #dominik code here
    return None


def test_api(seq, region, window, span, locws=None, locwe=None):

    ### Test with cutout after fold
    data_after = {'up': []}

    md = RNA.md()
    md.max_bp_span = span
    md.window_size = window
    # create new fold_compound object
    fc = RNA.fold_compound(str(seq), md, RNA.OPTION_WINDOW)
    # call prop window calculation
    fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_after)
    # fc.probs_window(region, RNA.PROBS_WINDOW_BPP, bpp_callback, data)

    if locws is not None and locwe is not None:  # If we only need a subset of the folded sequence
        data_after['up'] = [data_after['up'][x] for x in range(locws, locwe+1)]
        seq = seq[locws-1:locwe]


    ### Test with cutout before fold
    data_bef = {'up': []}

    if locws is not None and locwe is not None:  # If we only need a subset of the folded sequence
        seq = seq[locws-1:locwe]

    md = RNA.md()
    md.max_bp_span = span
    md.window_size = window
    # create new fold_compound object
    fc = RNA.fold_compound(str(seq), md, RNA.OPTION_WINDOW)
    # call prop window calculation
    fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_bef)
    # fc.probs_window(region, RNA.PROBS_WINDOW_BPP, bpp_callback, data)


def up_callback(v, v_size, i, maxsize, what, data):

    logid = scriptname+'.up_callback: '
    try:
        if what:
            data['up'].extend([v])
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))
