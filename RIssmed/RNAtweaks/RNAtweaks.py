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
from typing import Iterable, Tuple
from tempfile import NamedTemporaryFile, TemporaryDirectory
import subprocess
import numpy as np
import gzip
import math
from collections import defaultdict
import logging
import RNA

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


def isvalid(x=None):
    logid = scriptn+'.isvalid: '
    try:
        if x or x == 0:
            if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
                return False
            else:
                return True
        else:
            return False
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def isinvalid(x=None):
    logid = scriptn+'.isinvalid: '
    try:
        if x or x == 0:
            if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
                return True
            else:
                return False
        else:
            return True
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

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



def up_callback(v, v_size, i, maxsize, what, data):

    logid = scriptn + '.up_callback: '
    try:
        if what:
            data['up'].extend([v])
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


class PLFoldOutput:
    """Output wrapper for unpaired probability files of RNAplfold
           Attributes
          ----------
           text : str
              string representation of an RNAplfold punpaired output file.
              Missing the lines starting with # is supported

           """
    def __init__(self, text: str):
        self.text = self.__sanitize(text)
        self._numpy_array = None

    def __str__(self):
        return self.text

    def __eq__(self, other: PLFoldOutput):
        if self.__class__ != other.__class__:
            raise NotImplementedError
        return np.allclose(self.numpy_array, other.numpy_array, equal_nan=True, atol=0.0000001, rtol=0)

    @classmethod
    def from_file(cls, file_path: str):
        """creates PLfoldOutput from a punpaired file

            Parameters
           ----------
            file_path : str
               file location of the punpaired file
            Returns
           -------
           PLFoldOutput
               PLFoldOutput object
        """

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
        """creates PLfoldOutput from original rissmed np arrays

            Parameters
           ----------
            file_path : str
               file location of the rissmed numpy array
            Returns
           -------
           PLFoldOutput
               PLFoldOutput object
        """
        assert os.path.isfile(file_path), f"{file_path} is not a valid file path"
        ris_array = np.load(file_path)
        array = np.squeeze(ris_array)
        array_string = cls.__array_to_string(array)
        output = PLFoldOutput(array_string)
        output.__set_array(array)
        return output

    @classmethod
    def from_numpy(cls, array: np.ndarray):
        """creates PLFoldOutput from np arrays

            Parameters
           ----------
            array : np.ndarray
               array containing unpaired probabilities
            Returns
           -------
           PLFoldOutput
               PLFoldOutput object
        """
        array_string = cls.__array_to_string(array)
        output = PLFoldOutput(array_string)
        output.__set_array(array)
        return output

    @staticmethod
    def __sanitize(text: str):
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

    @staticmethod
    def __array_to_string(array):
        array = np.array(array, dtype=np.float).round(7)
        array_string = "\n".join(
            ['\t'.join([str(x + 1)] + [str(num) for num in array[x]]) for x in range(len(array))])
        return array_string

    def __set_array(self, array: np.ndarray):
        """sets numpy array and ensures data type is np.float"""
        array = np.array(array, dtype=np.float)
        self._numpy_array = array

    @property
    def numpy_array(self):
        """numpy array representation of lunpaired file"""
        if self._numpy_array is None:
            array = []
            for line in self.text.split("\n"):
                if not line.startswith("#") and not line.startswith(" #") and not line == "":
                    data = line.split("\t")[1:]
                    data = [float(x) if x != "NA" else np.nan for x in data]
                    array.append(data)
            array = np.array(array, dtype=np.float)
            self.__set_array(array)
        return self._numpy_array

    def localize(self, start: int, end: int):
        """trims the output (zero based)"""
        if self._numpy_array is not None:
            self._numpy_array = self._numpy_array[start:end]
        text_list = self.text.split("\n")[start+2:end+2]
        for x, item in enumerate(text_list):
            text_list[x] = "\t".join([str(x+1)] + item.split("\t")[1:])
        self.text = self.__sanitize("\n".join(text_list))

    def get_text(self, nan="NA", truncated=True) -> str:
        """get string of the unpaired probabiity file

            Parameters
           ----------
            nan : str, optional
               replaces not a number with this string (default is as original output 'NA')
            truncated : bool, optional
               choose whether the starting # lines should be included in the string
               (default is True)

            Returns
           -------
           str
               string representation of the unpaired probability output
        """
        out_string = self.text
        if truncated:
            out_string = "\n".join(out_string.split("\n")[2:])
        out_string = out_string.replace("NA", nan)
        return out_string

    def get_rissmed_np_array(self) -> np.ndarray:
        """get original RIssmed numpy array output

            Returns
           -------
           np.ndarray
               original rissmed numpy array with an additional empty axis
        """
        array = self.numpy_array
        array = np.array(array, dtype=np.float)
        array = np.expand_dims(array, axis=1)
        return array


def cmd_rnaplfold(sequence: str, window: int, span: int, region: int = 30, temperature: float = 37,
                  constraint: Iterable[Tuple[str, int, int]] = None) -> PLFoldOutput:
    """command line wrapper for RNAplfold

           Parameters
           ----------
            sequence : str
               string representation of the sequence either RNA or DNA
            window : int
               RNAplfold window option
            span: int
                RNAplfold span option
            region: int, optional
                RNAplfold region (u) option (default is 30)
            temperature: float
                RNAplfold temperature setting
            constraint: Iterable[Tuple[str, int, int]], optional
                Constraints as Tuple in format (paired(p)/unpaired(u), start, end) (default is None)
                !!Warning!! ZERO BASED !!Warning!!

           Returns
           -------
           PLFoldOutput
               PLFoldOutput object
           """
    with TemporaryDirectory() as tmp_dir, NamedTemporaryFile(mode="r+") as constraint_file:
        constraint_string = ""
        if constraint is not None:
            for entry in constraint:
                mode = entry[0]
                start = entry[1]
                end = entry[2]
                if mode == "paired" or mode == "p":
                    const = "F"
                elif mode == "unpaired" or mode == "u":
                    const = "P"
                else:
                    raise ValueError("Constraint wrongly formatted. Has to be ('paired(p)'/'unpaired(u)', start, end)")
                constraint_string += f"{const} {start+1} {0} {end - start}\n"
        constraint_file.write(constraint_string)
        constraint_file.seek(0)
        rnaplfold = subprocess.Popen(["RNAplfold", "-W", str(window), "-L", str(span),
                                      "--commands", constraint_file.name, "--auto-id", "-u", str(region), "-T",
                                      str(temperature)],
                                     stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                     cwd=tmp_dir)
        stdout, stderr = rnaplfold.communicate(sequence.encode("utf-8"))
        assert stderr == b"", f"call to RNApfold went wrong: \n {stderr.decode()}"
        file = os.path.join(tmp_dir, "sequence_0001_lunp")
        rnaplfold_output = PLFoldOutput.from_file(file)
        return rnaplfold_output


def api_rnaplfold(sequence: str, window: int, span: int, region: int = 30, temperature: float = 37,
                  constraint: Iterable[Tuple] = None) -> PLFoldOutput:
    """api wrapper for RNAplfold

           Parameters
           ----------
            sequence : str
               string representation of the sequence either RNA or DNA
            window : int
               RNAplfold window option
            span: int
                RNAplfold span option
            region: int, optional
                RNAplfold region (u) option (default is 30)
            temperature: float
                RNAplfold temperature setting
            constraint: Iterable[Tuple[str, int, int]], optional
                Constraints as Tuple in format (paired(p)/unpaired(u), start, end) (default is None)
                !!Warning!! ZERO BASED !!Warning!!

           Returns
           -------
           PLFoldOutput
               PLFoldOutput object
           """
    sequence = sequence.upper().replace("T", "U")
    data = {'up': []}
    md = RNA.md()
    md.max_bp_span = span
    md.window_size = window
    md.temperature = temperature

    # create new fold_compound object
    fc = RNA.fold_compound(str(sequence), md, RNA.OPTION_WINDOW)
    if constraint is not None:
        for entry in constraint:
            mode = entry[0]
            start = entry[1]
            end = entry[2]
            if mode == "paired" or mode == "p":
                fc = constrain_paired(fc, start, end)
            elif mode == "unpaired" or mode == "u":
                fc = constrain_unpaired(fc, start, end)
            else:
                raise ValueError("Constraint wrongly formatted. Has to be ('paired(p)'/'unpaired(u)', start, end)")

    # call prop window calculation
    fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data)
    array = np.array(data["up"]).squeeze()[:, 1:]
    pl_output = PLFoldOutput.from_numpy(array)
    return pl_output



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
