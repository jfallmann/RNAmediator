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
from __future__ import annotations
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
from lib.Collection import isinvalid

####################
# ViennaRNA helper
####################

try:
    log = logging.getLogger(__name__)  # use module name
    scriptn = __name__  # os.path.basename(inspect.stack()[-1].filename).replace('.py', '')
    log.debug('LOGGING IN RNAtweaks'+str(scriptn)+str(log)+str(log.handlers))
except Exception:
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
    except Exception:
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
    except Exception:
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

    except Exception:
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

    except Exception:
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
    except Exception:
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
    except Exception:
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

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def print_up(data=None, seqlength=None, region=None):
    logid = scriptn+'.print_up: '
    log.debug(logid+str(len(data))+' '+str(seqlength)+' '+str(region))
    try:
        if data:
            ups = ''
            if seqlength and seqlength != len(data):
                log.error(logid+'Lengths of sequence and array do not match: '+str(seqlength)+' and '+str(len(data)))
            for i in range(len(data)):
                if i >= len(data):
                    log.error(logid+'i larger than size of array')
                for x in range(1,region+1):
                    if x >= len(data[i]):
                        log.error(logid+'x larger than size of subarray')
                    if isinvalid(data[i][x]):
                        data[i][x] = np.nan
                    else:
                        data[i][x] = round(data[i][x],7)
                ups+=str(i+1)+"\t"+"\t".join(map(str,data[i][1:region+1]))+"\n"
            return ups
        else:
            log.error(logid+'No up data to print')
            return ups
    except Exception:
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
            return np.empty(region)
    except Exception:
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
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def printdiff(a, o=None):
    logid = scriptn+'.printdiff: '
    try:
        np.save(o, a)
    except Exception:
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
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def pl_to_array(name, ulim, fmt='npy'):
    logid = scriptn+'.pl_to_array: '
    try:
        log.debug('\t'.join([logid, name, str(ulim), fmt]))
        if fmt == 'txt':
            return np.array(np.loadtxt(name, usecols=ulim, unpack=True, delimiter='\t', encoding='bytes'))
        elif fmt == 'npy':
            #log.debug(np.load(name)[:,0][:,0])
            return np.array(np.load(name)[:,0][:,ulim-1])
    except Exception:
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
    except Exception:
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
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


class PLFoldOutput:
    def __init__(self, text: str):
        self.text = self._sanitize(text)

        self._array = None

    def __str__(self):
        return self.text

    def __eq__(self, other: PLFoldOutput):
        if self.__class__ != other.__class__:
            return False
        return np.array_equal(self.get_numpy_array(), other.get_numpy_array(), equal_nan=True)

    @staticmethod
    def _sanitize(text: str):
        list_text = text.split("\n")
        if list_text[0].startswith("1\t"):
            length = len(list_text[0].split("\t")) - 1
            part = [str(x + 1) for x in range(length)]
            list_text = ["#unpaired probabilities", " #$i\tl=" + "\t".join(part)] + list_text
            text = "\n".join(list_text)
        elif list_text[0].startswith("#") and list_text[1].startswith(" #"):
            pass
        else:
            raise ValueError("text seems not to be a valid plfold string")
        text = text.replace("nan", "NA")
        return text

    def get_numpy_array(self):
        if self._array is None:
            array = []
            for line in self.text.split("\n"):
                if not line.startswith("#") and not line.startswith(" #") and not line == "":
                    data = line.split("\t")[1:]
                    data = [float(x) if x != "NA" else np.nan for x in data]
                    array.append(data)
            array = np.array(array)
            self._array = array
        return self._array

    def set_array(self, array: np.ndarray):
        self._array = array

    def localize(self, start: int, end: int):
        if self._array is not None:
            self._array = self._array[start:end]
        text_list = self.text.split("\n")[start+2:end+2]
        for x, item in enumerate(text_list):
            text_list[x] = "\t".join([str(x+1)] + item.split("\t")[1:])
        self.text = self._sanitize("\n".join(text_list))

    @classmethod
    def from_file(cls, file_path: str):
        assert os.path.isfile(file_path), f"{file_path} is not a valid file path"
        if ".gz" in file_path:
            with gzip.open(file_path, "rt") as handle:
                file_data = handle.read()
        else:
            with open(file_path, "r") as handle:
                file_data = handle.read()
        return PLFoldOutput(file_data)

    @classmethod
    def from_rissmed_numpy_output(cls, file_path: str):
        assert os.path.isfile(file_path), f"{file_path} is not a valid file path"
        ris_array = np.load(file_path)
        array = np.squeeze(ris_array)
        array_string = cls._array_to_string(array)
        output = PLFoldOutput(array_string)
        output.set_array(array)
        return output

    @classmethod
    def from_numpy(cls, array: np.ndarray):
        array_string = cls._array_to_string(array)
        output = PLFoldOutput(array_string)
        output.set_array(array)
        return output

    @staticmethod
    def _array_to_string(array):
        array = np.array(array, dtype=np.float).round(7)
        array_string = "\n".join(
            ['\t'.join([str(x + 1)] + [str(num) for num in array[x]]) for x in range(len(array))])
        return array_string

    def get_text(self, nan="NA", truncated=True):
        out_string = self.text
        if truncated:
            out_string = "\n".join(out_string.split("\n")[2:])
        out_string = out_string.replace("NA", nan)
        return out_string

    def get_rissmed_np_array(self):
        array = self.get_numpy_array()
        array = np.array(array, dtype=np.float)
        array = np.expand_dims(array, axis=1)
        return array





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
