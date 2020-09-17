# RNAtweaks.py ---
#
# Filename: RNAtweaks.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Fri Aug 21 10:23:40 2020 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Sep  1 11:09:40 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 16
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
import traceback as tb
import numpy as np
import gzip
import math
from collections import defaultdict
import logging
# own
from lib.Collection import *

####################
# ViennaRNA helper
####################

try:
    log = logging.getLogger(__name__)  # use module name
    scriptn = os.path.basename(inspect.stack()[-1].filename).replace('.py', '')
    log.debug('LOGGING IN RNAtweaks'+str(scriptn)+str(log)+str(log.handlers))
except Exception as err:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    print(''.join(tbe.format()),file=sys.stderr)

### Calculate nrg/prob/bppm/ddg

def calc_gibbs(fc):
    logid = scriptn+'.calc_gibbs: '
    try:
        return fc.pf()[1]
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def get_bppm(tmp, start, end):
    logid = scriptn+'.get_bppm: '
    bppm = []
    try:
        if start < 0 or end > len(tmp):
            log.warning(logid+'start of constraint '+str(start)+' end of constraint '+str(end)+' while length of bpp matrix '+str(len(tmp))+'! Skipping!')
            return None

        for item in tmp:
            for i in range(int(start),int(end)+1):
                if item[i] > 0.0:
                    bppm.append(str.join('\t',[str(tmp.index(item)), str(i), str(item[i])]))
        return bppm
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def get_ddg(file):
    logid = scriptn+'.get_ddg: '
    try:
        ret = defaultdict()
        if (isinstance(file, str) and os.path.isfile(file)):
            if '.gz' in file :
                res = gzip.open(file,'rt')
            else:
                res = open(file,'rt')

            for line in res:
                log.debug(logid+line)
                if 'Condition' in line[0:15]:
                    continue
                else:
                    cond, gibbs, dg, nrg, cons = line.rstrip().split('\t')
                    if not str(cons) in ret:
                        ret[str(cons)] = defaultdict()
                    ret[str(cons)][cond] = float(dg)
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def calc_ddg(ddgs):
    logid = scriptn+'.calc_ddg: '

    try:
        log.debug(logid+str(ddgs))
        ddg = ddgs['constraint_unpaired']+ddgs['secondconstraint_unpaired']-ddgs['bothconstraint_unpaired']-ddgs['unconstraint']
        """Yi-Hsuan Lin, Ralf Bundschuh, RNA structure generates natural cooperativity between single-stranded RNA binding proteins targeting 5' and 3'UTRs, Nucleic Acids Research, Volume 43, Issue 2, 30 January 2015, Pages 1160-1169, https://doi.org/10.1093/nar/gku1320"""

        return ddg

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def calc_bpp(bppm):
    logid = scriptn+'.calc_bpp: '
    bpp = 0.0;
    try:
        for entry in bppm:
            base, mate, prob = map(float,entry.split('\t'))
            bpp += prob
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

    return bpp

def calc_nrg(bpp):
    logid = scriptn+'.calc_nrg: '
    #set kT for nrg2prob and vice versa calcs
    kT = 0.61632077549999997
    nrg = 0.0;
    try:
        if bpp > 0.0:
            nrg = -1 * kT * math.log(bpp)
        return nrg
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def print_region_up(data, seqlength=None, region=None):
    logid = scriptn+'.print_region_up: '
    try:
        if data:
            ups=''
            x = int(region)
            for i in range(int(seqlength)):
                if isinvalid(data[i][x]):
                    data[i][x] = np.nan
                else:
                    data[i][x] = round(data[i][x],7)
                ups+=str(i+1)+"\t"+str(data[i][x])+"\n"
            return ups
        else:
            log.error(logid+'No up data to print')
            return ups

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def print_up(data=None, seqlength=None, region=None):
    logid = scriptn+'.print_up: '
    try:
        if data:
            ups=''
            for i in range(int(seqlength)):
                for x in range(1,region+1):
                    if isinvalid(data[i][x]):
                        data[i][x] = np.nan
                    else:
                        data[i][x] = round(data[i][x],7)
                ups+=str(i+1)+"\t"+"\t".join(map(str,data[i][1:region+1]))+"\n"
            return ups
        else:
            log.error(logid+'No up data to print')
            return ups
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def up_to_array(data=None, region=None, seqlength=None):
    logid = scriptn+'.up_to_array: '
    try:
        if data:
            entries=[]
            if not seqlength:
                seqlength = len(data)
            if not region:
                region = slice(1,len(data[0]))
            for i in range(seqlength):
                entries.append([])
                for e in range(len(data[i])):
                    if isinvalid(data[i][e]):
                        data[i][e] = np.nan
                    else:
                        data[i][e] = round(data[i][e],8)
                entries[i].append(data[i][region])
            return np.array(entries)
        else:
            log.error(logid+'No up data to print')
            return np.array()
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format())+'\t'+str(entries)+'\t'+str(region)+'\t'+str(seqlength))

def npprint(a, o=None):#, format_string ='{0:.2f}'):
    logid = scriptn+'.npprint: '
    try:
        out = ''
        it = np.nditer(a, flags=['f_index'])
        while not it.finished:
            out += "%d\t%0.7f" % (it.index+1,it[0])+"\n"
            it.iternext()
        if o:
            o.write(bytes(out,encoding='UTF-8'))
        else:
            print(out)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def printdiff(a, o=None):
    logid = scriptn+'.printdiff: '
    try:
        np.save(o, a)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def read_precalc_plfold(data, name, seq):
    logid = scriptn+'.read_precalc_plfold: '
    try:
        for i in range(len(seq)):
            data.append([])
            data[i] = []
        with gzip.open(name,'rt') as o:
            for line in o:
                cells = line.rstrip().split('\t')
                data[int(cells[0])-1].append([])
                data[int(cells[0])-1][0] = None
                for a in range(1,len(cells)):
                    data[int(cells[0])-1].append([])
                    data[int(cells[0])-1][a] = float(cells[a])
        return data
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def pl_to_array(name, ulim, fmt='npy'):
    logid = scriptn+'.pl_to_array: '
    try:
        log.debug('\t'.join([logid,name]))
        if fmt == 'txt':
            return np.array(np.loadtxt(name, usecols=ulim, unpack=True, delimiter='\t', encoding='bytes'))
        elif fmt == 'npy':
            return np.array(np.load(name)[:,ulim-1])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+' '+name+': '.join(tbe.format()))

###Constraints

def constrain_paired(fc, start, end):
    logid = scriptn+'.constrain_paired: '
    try:
        for x in range(start+1, end+1):
            fc.hc_add_bp_nonspecific(x,0) #0 means without direction  ( $ d < 0 $: pairs upstream, $ d > 0 $: pairs downstream, $ d == 0 $: no direction)
        return fc
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

def constrain_unpaired(fc, start, end):
    logid = scriptn+'.constrain_unpaired: '
    try:
        for x in range(start+1, end+1):
            fc.hc_add_up(x)
        return fc
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


#def bpp_callback(v, v_size, i, maxsize, what, data):
#   if what & RNA.PROBS_WINDOW_BPP:
#       data['bpp'].extend([{'i': i, 'j': j, 'p': p} for j, p in enumerate(v) if (p is not None)])# and (p >= 0.01)])
#
#def up_callback(v, v_size, i, maxsize, what, data):
#   if what & RNA.PROBS_WINDOW_UP:
#       #    data['up'].extend([{ 'i': i, 'up': v}])
#       data['up'].extend([v])


#
# RNAtweaks.py ends here
