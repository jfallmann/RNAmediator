from __future__ import annotations

import gzip
import logging
import math
import os
import subprocess
import sys
import traceback as tb
from collections import defaultdict
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Iterable, Tuple

# set path for VRNA lib if necessary
# not supported anymore needs to be added in the RNAtweaks file
# if vrna:
#     sys.path = [vrna] + sys.path
#     global RNA
#     RNA = importTweaks.import_module('RNA')
#     globals().update(
#         {n: getattr(RNA, n) for n in RNA.__all__}
#         if hasattr(RNA, '__all__')
#         else {k: v for (k, v) in RNA.__dict__.items() if not k.startswith('_')}
#        )

import RNA
import numpy as np

####################
# ViennaRNA helper
####################

try:
    log = logging.getLogger(__name__)  # use module name
    scriptn = (
        __name__  # os.path.basename(inspect.stack()[-1].filename).replace('.py', '')
    )
    log.debug("LOGGING IN RNAtweaks" + str(scriptn) + str(log) + str(log.handlers))
except Exception:
    EXC_TYPE, EXC_VALUE, EXC_TB = sys.exc_info()
    TBE = tb.TracebackException(
        EXC_TYPE,
        EXC_VALUE,
        EXC_TB,
    )
    print("".join(TBE.format()), file=sys.stderr)


def _isvalid(x=None):
    """Checks if x is a valid value

    Parameters
    ----------
     x: any
        The value to check

    Returns
    -------
    bool
        True/False
    """

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


def _isinvalid(x=None):
    """Checks if x is invalid

    Parameters
    ----------
     x: any
        The value to check

    Returns
    -------
    bool
        True/False
    """

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


# Calculate nrg/prob/bppm/ddg


def _calc_gibbs(fc):
    """Calcuates Gibbs Free Energy for sequence

    Parameters
    ----------
     fc: object
        RNA.fold_compound object of ViennaRNA API with or without added constraints

    Returns
    -------
    fc.pf()[1]: float
        Partition Function free energy value
    """

    logid = scriptn + ".calc_gibbs: "
    try:
        return fc.pf()[1]
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _get_bppm(matrix, start, end):
    """Extracts base pair probability matrix for sequence from start to end

    Parameters
    ----------
    matrix: object
        RNA.fold_compound, BPPM of structure calculated by ViennaRNA
    start: int
    end: int

    Returns
    -------
    bppm: object
        RNA.fold_compound subset of original bppm
    """

    logid = scriptn + ".get_bppm: "
    bppm = []
    try:
        if start < 0 or end > len(tmp):
            log.warning(
                f"{logid} start of constraint: {start} end of constraint: {end} "
                f"while length of bpp matrix {len(tmp)} ! Skipping!"
            )
            return None

        for item in tmp:
            for i in range(int(start), int(end) + 1):
                if item[i] > 0.0:
                    bppm.append(
                        str.join("\t", [str(tmp.index(item)), str(i), str(item[i])])
                    )
        return bppm
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _get_ddg(file):
    """Derives Delta-Delta Gibbs Free Energy for sequence from output file

    Parameters
    ----------
     file: str
        filename

    Returns
    -------
    ret: defaultdict[str, defaultdict[str, float]]
        Partition Function free energy value
    """

    logid = scriptn + ".get_ddg: "
    try:
        ret = defaultdict()
        if isinstance(file, str) and os.path.isfile(file):
            if ".gz" in file:
                res = gzip.open(file, "rt")
            else:
                res = open(file, "rt")

            for line in res:
                log.debug(logid + line)
                if "Condition" in line[0:15]:
                    continue
                else:
                    cond, gibbs, dg, nrg, cons = line.rstrip().split("\t")
                    if not str(cons) in ret:
                        ret[str(cons)] = defaultdict()
                    ret[str(cons)][cond] = float(dg)
        return ret

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _calc_ddg(ddgs):
    """Calculates Delta-Delta-Gibbs Free Energy

    Args:
    ddgs: dict
        Dictionary containing relevant information for ddg calculation, produced by _get_ddgs

    Returns:
    ddg: float
        Delta-Delta-Gibbs Free Energy after binding of two antagonistic/cooperative binding factors

    Reference:
        Yi-Hsuan Lin, Ralf Bundschuh, RNA structure generates natural cooperativity between single-stranded RNA binding proteins targeting 5' and 3'UTRs, Nucleic Acids Research, Volume 43, Issue 2, 30 January 2015, Pages 1160-1169, https://doi.org/10.1093/nar/gku1320
    """

    logid = scriptn + ".calc_ddg: "

    try:
        log.debug(logid + str(ddgs))
        cons_up = ddgs["constraint_unpaired"]
        sec_cons_up = ddgs["secondconstraint_unpaired"]
        both_cons_up = ddgs["bothconstraint_unpaired"]
        uncons = ddgs["unconstraint"]
        ddg = cons_up + sec_cons_up - both_cons_up - uncons

        return ddg

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _calc_bpp(bppm):
    """Calculate base-pair probability

    Parameters
    ----------
    bppm : list
        Base-pair probability matrix derived from ViennaRNA-API

    Returns
    -------
    bpp: float
        base-pair probability
    """

    logid = scriptn + ".calc_bpp: "
    bpp = 0.0
    try:
        for entry in bppm:
            base, mate, prob = map(float, entry.split("\t"))
            bpp += prob
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))

    return bpp


def _calc_nrg(bpp):
    """Calculate pseudo opening energy from base-pair probability

    Parameters
    ----------
    bpp : float
        Base-pair probabilit, derived from e.g. _calc_bpp

    Returns
    -------
    nrg: float
        Pseudo opening energy
    """

    logid = scriptn + ".calc_nrg: "
    # set k_t for nrg2prob and vice versa calcs
    k_t = 0.61632077549999997
    nrg = 0.0
    try:
        if bpp > 0.0:
            nrg = -1 * k_t * math.log(bpp)
        return nrg
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _print_region_up(data, seqlength=None, region=None):
    """Print pairing probability of region

    Parameters
    ----------
    data : dict[int, dict[int]]
        Output of ViennaRNA fold_compound folding
    seqlength : int, optional
        Length of sequence to print, by default None
    region : int, optional
        Position in data dict, corressponds to -u parameter of RNAplfold, by default None

    Returns
    -------
    ups: str
        String of propabilities of being unpaired for given region

    Raises
    ------
    NameError
        If no data dict was found
    """

    logid = scriptn + ".print_region_up: "
    try:
        if data:
            ups = ""
            x = int(region)
            for i in range(int(seqlength)):
                if _isinvalid(data[i][x]):
                    data[i][x] = np.nan
                else:
                    data[i][x] = round(data[i][x], 7)
                ups += str(i + 1) + "\t" + str(data[i][x]) + "\n"
            return ups
        else:
            log.error(logid + "No up data to print")
            raise NameError("name 'ups' is not defined")
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _print_up(data=None, seqlength=None, region=None):
    """Print pairing probability

    Parameters
    ----------
    data : dict[int, dict[int]]
        Output of ViennaRNA fold_compound folding
    seqlength : int, optional
        Length of sequence to print, by default None
    region : int, optional
        Position in data dict, corressponds to -u parameter of RNAplfold, by default None

    Returns
    -------
    ups: str
        String of propabilities of being unpaired for sequence

    Raises
    ------
    NameError
        If no data dict was found
    """

    logid = scriptn + ".print_up: "
    log.debug(logid + str(len(data)) + " " + str(seqlength) + " " + str(region))
    try:
        if data:
            ups = ""
            if seqlength and seqlength != len(data):
                log.error(
                    logid
                    + "Lengths of sequence and array do not match: "
                    + str(seqlength)
                    + " and "
                    + str(len(data))
                )
            for i in range(len(data)):
                if i >= len(data):
                    log.error(logid + "i larger than size of array")
                for x in range(1, region + 1):
                    if x >= len(data[i]):
                        log.error(logid + "x larger than size of subarray")
                    if _isinvalid(data[i][x]):
                        data[i][x] = np.nan
                    else:
                        data[i][x] = round(data[i][x], 7)
                ups += (
                    str(i + 1)
                    + "\t"
                    + "\t".join(map(str, data[i][1 : region + 1]))
                    + "\n"
                )
            return ups
        else:
            log.error(logid + "No up data to print")
            raise NameError("name 'ups' is not defined")
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _up_to_array(data=None, region=None, seqlength=None):
    """Converts output of ViennaRNA folc_compound.probs_window to numpy array

    Parameters
    ----------
    data : dict[int,dict[int]], optional
        Dictionary holding return value of robs_window, by default None
    region : int, optional
        Region of interest, corresponding to -u parameter or RNAPLfold, by default None
    seqlength : int, optional
        Length of sequence that was folded, by default None

    Returns
    -------
    numpy.array
        Numpy array holding probabilities of being unpaired for region of interest
    """

    logid = scriptn + ".up_to_array: "
    entries = []
    try:
        if data:
            if not seqlength:
                seqlength = len(data)
            if not region:
                region = slice(1, len(data[0]))
            for i in range(seqlength):
                entries.append([])
                for e in range(len(data[i])):
                    if _isinvalid(data[i][e]):
                        data[i][e] = np.nan
                    else:
                        data[i][e] = round(data[i][e], 8)
                entries[i].append(data[i][region])
            return np.array(entries)
        else:
            log.error(logid + "No up data to print")
            return np.empty(region)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(
            logid
            + "".join(tbe.format())
            + "\t"
            + str(entries)
            + "\t"
            + str(region)
            + "\t"
            + str(seqlength)
        )


def _npprint(a, o=None):  # format_string ='{0:.2f}'):
    """Pretty print of numpy array"""

    logid = scriptn + ".npprint: "
    try:
        out = ""
        it = np.nditer(a, flags=["f_index"])
        while not it.finished:
            out += "%d\t%0.7f" % (it.index + 1, it[0]) + "\n"
            it.iternext()
        if o:
            o.write(bytes(out, encoding="UTF-8"))
        else:
            print(out)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def printdiff(a, o=None):
    """Print difference between numpy arrays

    Parameters
    ----------
    a : array
        Array to print
    o : str, optional
        where to save the diff to, by default None
    """

    logid = scriptn + ".printdiff: "
    try:
        np.save(o, a)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _read_precalc_plfold(data, name, seq):
    """Reads in plfold output of previous run

    Parameters
    ----------
    data : list
        Where to store the data to
    name : str
        Name of file to read from
    seq : str
        Sequence that was folded

    Returns
    -------
    data: list[int, list[int]]
        List containing precalculated results
    """

    logid = scriptn + ".read_precalc_plfold: "
    try:
        for i in range(len(seq)):
            data.append([])
            data[i] = []
        with gzip.open(name, "rt") as o:
            for line in o:
                cells = line.rstrip().split("\t")
                data[int(cells[0]) - 1].append([])
                data[int(cells[0]) - 1][0] = None
                for a in range(1, len(cells)):
                    data[int(cells[0]) - 1].append([])
                    data[int(cells[0]) - 1][a] = float(cells[a])
        return data
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _pl_to_array(name, ulim, fmt="npy"):
    """Tranforms RNAPLFold output to numpy array

    Parameters
    ----------
    name : str
        Name of file to read in
    ulim : int
        Region of interest, corresponds to -u parameter of RNAPLfold
    fmt : str, optional
        Format of input file, by default "npy"

    Returns
    -------
    numpy.array
        Numpy array containing data from read in file
    """

    logid = scriptn + ".pl_to_array: "
    try:
        log.debug("\t".join([logid, name, str(ulim), fmt]))
        if fmt == "txt":
            return np.array(
                np.loadtxt(
                    name, usecols=ulim, unpack=True, delimiter="\t", encoding="bytes"
                )
            )
        elif fmt == "npy":
            # log.debug(np.load(name)[:,0][:,0])
            return np.array(np.load(name)[:, 0][:, ulim - 1])
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + " " + name + ": ".join(tbe.format()))


def get_location(entry):
    """Returns location of folding window in genomic or local coordinates

    Parameters
    ----------
    entry : str
        String containing start, end and strand of folding window

    Returns
    -------
    ret: list(list(int,int,str))
        List of start, end and strand for folding window
    """

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


# Constraints


def _constrain_paired(fc, start, end):
    """Adds hard constraint paired to region

    Parameters
    ----------
    fc : RNA.fold_compound object
        ViennaRNA fold_compound object
    start : int
        Start of constraint
    end : int
        End of constraint

    Returns
    -------
    fc: RNA.fold_compound object
        ViennaRNA fold_compound object with constraint
    """

    logid = scriptn + ".constrain_paired: "
    try:
        for x in range(start + 1, end + 1):
            # 0 means without direction
            # ( $ d < 0 $: pairs upstream, $ d > 0 $: pairs downstream, $ d == 0 $: no direction)
            fc.hc_add_bp_nonspecific(x, 0)
        return fc
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _constrain_unpaired(fc, start, end):
    """Adds hard constraint unpaired to region

    Parameters
    ----------
    fc : RNA.fold_compound object
        ViennaRNA fold_compound object
    start : int
        Start of constraint
    end : int
        End of constraint

    Returns
    -------
    fc: RNA.fold_compound object
        ViennaRNA fold_compound object with constraint
    """

    logid = scriptn + ".constrain_unpaired: "
    try:
        for x in range(start + 1, end + 1):
            fc.hc_add_up(x)
        return fc
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def bpp_callback(v, v_size, i, maxsize, what, data):
    """Callback function for ViennaRNA API base-pair probability calculation

    Parameters
    ----------
    v : ViennaRNA fold_compound object
        The fold compound to use
    v_size : int
        size of the compound
    i : int
        From
    maxsize : int
        To
    what : str
        Which type of fold_compound
    data : dict[list(int, int, int)]
        Datastructure to fill and return
    """

    logid = scriptn + ".bpp_callback: "
    try:
        if what:
            data["bpp"].extend(
                [{"i": i, "j": j, "p": p} for j, p in enumerate(v) if (p is not None)]
            )  # and (p >= 0.01)])
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def _up_callback(v, v_size, i, maxsize, what, data):
    """Callback function for ViennaRNA API probability of being unpaired calculation

    Parameters
    ----------
    v : ViennaRNA fold_compound object
        The fold compound to use
    v_size : int
        size of the compound
    i : int
        From
    maxsize : int
        To
    what : str
        Which type of fold_compound
    data : dict[list(int, int, int)]
        Datastructure to fill and return
    """

    logid = scriptn + ".up_callback: "
    try:
        if what:
            data["up"].extend([v])
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


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
        return np.allclose(
            self.numpy_array, other.numpy_array, equal_nan=True, atol=0.0000001, rtol=0
        )

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
            list_text = [
                "#unpaired probabilities",
                " #$i\tl=" + "\t".join(part),
            ] + list_text
            text = "\n".join(list_text)
        elif list_text[0].startswith("#") and list_text[1].startswith(" #"):
            pass
        else:
            raise ValueError("text seems not to be a valid plfold string")
        text = text.replace("nan", "NA")
        return text

    @staticmethod
    def __array_to_string(array):
        array = np.array(array, dtype=float).round(7)
        array_string = "\n".join(
            [
                "\t".join([str(x + 1)] + [str(num) for num in array[x]])
                for x in range(len(array))
            ]
        )
        return array_string

    def __set_array(self, array: np.ndarray):
        """sets numpy array and ensures data type is float"""
        array = np.array(array, dtype=float)
        self._numpy_array = array

    @property
    def numpy_array(self):
        """numpy array representation of lunpaired file"""
        if self._numpy_array is None:
            array = []
            for line in self.text.split("\n"):
                if (
                    not line.startswith("#")
                    and not line.startswith(" #")
                    and not line == ""
                ):
                    data = line.split("\t")[1:]
                    data = [float(x) if x != "NA" else np.nan for x in data]
                    array.append(data)
            array = np.array(array, dtype=float)
            self.__set_array(array)
        return self._numpy_array

    def localize(self, start: int, end: int):
        """trims the output (zero based)"""
        if self._numpy_array is not None:
            self._numpy_array = self._numpy_array[start:end]
        text_list = self.text.split("\n")[start + 2 : end + 2]
        for x, item in enumerate(text_list):
            text_list[x] = "\t".join([str(x + 1)] + item.split("\t")[1:])
        self.text = self.__sanitize("\n".join(text_list))

    def get_text(self, nan="NA", truncated=True) -> str:
        """get string of the unpaired probability file

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
        array = np.array(array, dtype=float)
        array = np.expand_dims(array, axis=1)
        return array


class FoldOutput(defaultdict):
    """Output wrapper for ouput files of RNAfold

     Attributes
    ----------
     text : str
        string representation of an RNAfold output file.
     anno : str, optional
        Annotation for file
    """

    def __init__(self, text: str, condition: str = 'unconstraint', gibbs: float = None, constraint: str = None, bpp: float = None, bppm: np.array = None, nrg: float = None, ddg: float = None, dnrg: float = None):
        """Output wrapper for RNAfold

        Parameters
        ----------
        text : str
            string representation of an RNAfold output file
        condition : str, optional
            Condition used for folding, by default 'raw'
        gibbs : float, optional
            Gibbs Free Energy, by default None
        constraint : str, optional
            Constraint applied for folding, by default None
        bpp : float, optional
            Base-Pairing Probability, by default None
        bppm : np.array, optional
            BasePairProbability matrix, by default None
        nrg : float, optional
            Folding Pseudo-Energy, by default None
        ddg : float, optional
            Delta-delta-Gibbs to unconstraint, by default None
        dnrg : float, optional
            Delta-Folding Pseudo-Energy, by default None
        """

        assert any(condition == x for x in ["unconstraint", "constraint", "pairedconstraint", "constraint_unpaired", "constraint_paired", "secondconstraint_unpaired", "secondconstraint_paired", "bothconstraint_unpaired", "bothconstraint_paired"])

        self[condition] = defaultdict
        self[condition]['text'] = text
        self[condition]['gibbs'] = gibbs
        self[condition]['constraint'] = None
        self[condition]['bppm'] = np.array()
        self[condition]['bpp'] = None
        self[condition]['nrg'] = None
        self[condition]['ddg'] = None
        self[condition]['dnrg'] = None


    def __eq__(self, other: FoldOutput):
        if self.__class__ != other.__class__:
            raise NotImplementedError
        if np.assert_almost_equal(self._bppm, other._bppm) and np.allclose(
            self.bppm, other.bppm, equal_nan=True, atol=0.0000001, rtol=0
        ):
            return True
        else:
            return False

    
    @classmethod
    def from_RNAfold_file(cls, file_path: str, condition: str = 'raw', constraint: str = None):
        """creates FoldOutput from a file

         Parameters
        ----------
         file_path : str
            file location of the output file
         condition : str, optional
            which folding condition, e.g. raw, constraint, constraint-paired
         constraint : str, optional
            which constraint was applied when folding

         Returns
        -------
        FoldOutput
            FoldOutput object
        """

        assert os.path.isfile(file_path), f"{file_path} is not a valid file path"
        if ".gz" in file_path:
            with gzip.open(file_path, "rt") as handle:
                file_data = handle.readlines()
                gibbs = file_data[2].split(' ')[1].strip("()")
        else:
            with open(file_path, "r") as handle:
                file_data = handle.readlines()
                gibbs = file_data[2].split(' ')[1].strip("()")

        return FoldOutput(file_data, condition, gibbs)


    @classmethod
    def from_rissmed_fold_compound(cls, fc, condition: str = 'unconstraint', start: int, end: int):
        """creates FoldOutput from original rissmed fold compound

         Parameters
        ----------
         fc : fold_compound object
            ViennaRNA fold_compound object
         condition: str
            Condition for folding
         start: int
            Start of folding window
         end: int
            End of folding window

         Returns
        -------
        FoldOutput
            FoldOutput object
        """

        assert any(condition == x for x in ["unconstraint", "constraint", "pairedconstraint", "constraint_unpaired", "constraint_paired", "secondconstraint_unpaired", "secondconstraint_paired", "bothconstraint_unpaired", "bothconstraint_paired"])

        cls[condition]['gibbs'] = calc_gibbs(fc)
        cls[condition]['bppm'] = np.array(get_bppm(fc.bpp(), start, end))
        cls[condition]['bpp'] = calc_bpp(cls[condition]['bppm'])
        cls[condition]['nrg'] = calc_nrg(bpp)

        if cls.get('unconstraint') and condition != 'unconstraint':
            cls[condition]['ddg'] = cls['unconstraint']['gibbs'] - cls[condition]['gibbs']
            cls[condition]['dnrg'] = cls['unconstraint']['nrg'] - cls[condition]['nrg']
        else:
            cls[condition]['ddg'] = None
            cls[condition]['dnrg'] = None

        return cls


    @staticmethod
    def __sanitize(text: str):
        list_text = text.split("\n")
        if list_text[0].startswith("1\t"):
            length = len(list_text[0].split("\t")) - 1
            part = [str(x + 1) for x in range(length)]
            list_text = [
                "#unpaired probabilities",
                " #$i\tl=" + "\t".join(part),
            ] + list_text
            text = "\n".join(list_text)
        elif list_text[0].startswith("#") and list_text[1].startswith(" #"):
            pass
        else:
            raise ValueError("text seems not to be a valid plfold string")
        text = text.replace("nan", "NA")
        return text


    @staticmethod
    def __array_to_string(array):
        array = np.array(array, dtype=float).round(7)
        array_string = "\n".join(
            [
                "\t".join([str(x + 1)] + [str(num) for num in array[x]])
                for x in range(len(array))
            ]
        )
        return array_string


    def __set_array(self, condition: str = 'unconstraint', array: np.ndarray):
        """sets numpy array and ensures data type is float"""
         assert any(condition == x for x in ["unconstraint", "constraint", "pairedconstraint", "constraint_unpaired", "constraint_paired", "secondconstraint_unpaired", "secondconstraint_paired", "bothconstraint_unpaired", "bothconstraint_paired"])

        array = np.array(array, dtype=float)
        self[condition]["bppm"] = array


    def get_text(self):
        """return formatted string output

        Returns
        -------
        formatted_string
            Contents of FoldOutput as formatted string
        """
        
        formatted_string = list()
        formatted_string.append(str.join("\t", [
                                    "Condition",
                                    "FreeNRG(gibbs)",
                                    "deltaG",
                                    "OpeningNRG",
                                    "Constraint",
                                ],
                            )
                            + "\n",)) 
        for condition in self.keys():
            gibbs, ddg, nrg, const = [self[condition][x] for x in ['gibbs', 'ddg'. 'nrg', 'constraint']]
            formatted_string.append(str.join("\t", [condition, str(gibbs), str(ddg), str(nrg), str(const)]) + "\n")
        return formatted_string


def cmd_rnaplfold(
    sequence: str,
    window: int,
    span: int,
    region: int = 30,
    temperature: float = 37,
    constraint: Iterable[Tuple[str, int, int]] = None,
) -> PLFoldOutput:
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
    with TemporaryDirectory() as tmp_dir, NamedTemporaryFile(
        mode="r+"
    ) as constraint_file:
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
                    raise ValueError(
                        "Constraint wrongly formatted. Has to be ('paired(p)'/'unpaired(u)', start, end)"
                    )
                constraint_string += f"{const} {start+1} {0} {end - start}\n"
        constraint_file.write(constraint_string)
        constraint_file.seek(0)
        rnaplfold = subprocess.Popen(
            [
                "RNAplfold",
                "-W",
                str(window),
                "-L",
                str(span),
                "--commands",
                constraint_file.name,
                "--auto-id",
                "-u",
                str(region),
                "-T",
                str(temperature),
            ],
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            cwd=tmp_dir,
        )
        stdout, stderr = rnaplfold.communicate(sequence.encode("utf-8"))
        assert stderr == b"", f"call to RNApfold went wrong: \n {stderr.decode()}"
        file = os.path.join(tmp_dir, "sequence_0001_lunp")
        rnaplfold_output = PLFoldOutput.from_file(file)
        return rnaplfold_output


def api_rnaplfold(
    sequence: str,
    window: int,
    span: int,
    region: int = 5,
    temperature: float = 37,
    constraint: Iterable[Tuple] = None,
) -> PLFoldOutput:
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
    data = {"up": []}
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
                fc = _constrain_paired(fc, start, end)
            elif mode == "unpaired" or mode == "u":
                fc = _constrain_unpaired(fc, start, end)
            else:
                raise ValueError(
                    "Constraint wrongly formatted. Has to be ('paired(p)'/'unpaired(u)', start, end)"
                )

    # call prop window calculation
    fc.probs_window(region, RNA.PROBS_WINDOW_UP, _up_callback, data)
    array = np.array(data["up"]).squeeze()[:, 1:]
    pl_output = PLFoldOutput.from_numpy(array)
    return pl_output


def api_rnafold(
    sequence: str,
    window: int,
    span: int,
    temperature: float = 37,
    constraint: Iterable[Tuple] = None,
    FoldOut: FoldOutput object = None
) -> FoldOutput:
    """api wrapper for RNAfold

    Parameters
    ----------
     sequence : str
        string representation of the sequence either RNA or DNA
     window : int
        RNAfold window option
     span: int
         RNAfold basepair span option
     temperature: float
         RNAplfold temperature setting
     constraint: Iterable[Tuple[str, int, int, int, int]], optional
         Constraints as Tuple in format (paired(p)/unpaired(u), start, end) (default is None)
         !!Warning!! ZERO BASED !!Warning!!
     FoldOut: FoldOutput object, optional
        existing FoldOutput object to append use

    Returns
    -------
    FoldOutput
        FoldOutput object
    """

    # data
    data = {"seq": seqtofold, "stru": []}

    # set model details
    md = RNA.md()
    md.max_bp_span = span
    md.temperature = temperature

    # create new fold_compound object
    fc = RNA.fold_compound(data["seq"], md)

    if constraint is not None:
        for entry in constraint:
            mode = entry[0]
            s, e, os, oe = [None, None, None, None]
            start = entry[1]
            end = entry[2]
            fstart = entry[3] if len(entry > 2) else None
            fend = entry[4] if len(entry > 3) else None
            if any(mode == x for x in ["paired", "p", "constraint_paired"]): 
                fc = _constrain_paired(fc, start, end)
            elif mode == "secondconstraint_paired":
                fc = _constrain_paired(fc, fstart, fend)
            elif mode == "bothconstraint_paired":
                fc = _constrain_paired(fc, start, end)
                fc = _constrain_paired(fc, fstart, fend)
            elif any(mode == x for x in ["unpaired", "u", "constraint_unpaired"]): 
                fc = _constrain_unpaired(fc, start, end)
            elif mode == "secondconstraint_unpaired":
                fc = _constrain_unpaired(fc, fstart, fend)
            elif mode == "bothconstraint_unpaired":
                fc = _constrain_unpaired(fc, start, end)
                fc = _constrain_unpaired(fc, fstart, fend)
            else:
                raise ValueError(
                    "Constraint wrongly formatted. Has to be ('paired(p)'/'unpaired(u)', start, end)"
                )
    if not FoldOut:
        FoldOut = FoldOutput('none', mode)
    # call pf and prop calculation
    FoldOut[mode]['gibbs'] = fc.pf()[1]
    FoldOut[mode]['bppm'] = get_bppm(fc.bpp(), cstart, cend)

    if bppm is None:
        log.error(logid + "Empty bpp matrix returned, stopping here!")

    return FoldOut


# def bpp_callback(v, v_size, i, maxsize, what, data):
#   if what & RNA.PROBS_WINDOW_BPP:
#       data['bpp'].extend([{'i': i, 'j': j, 'p': p} for j, p in enumerate(v) if (p is not None)])# and (p >= 0.01)])
#
# def up_callback(v, v_size, i, maxsize, what, data):
#   if what & RNA.PROBS_WINDOW_UP:
#       #    data['up'].extend([{ 'i': i, 'up': v}])
#       data['up'].extend([v])


#
# Tweaks.py ends here
