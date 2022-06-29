#!/usr/bin/env python3
### GenerateBigWig.py ---
##
## Filename: GenerateBigWig.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Sep  6 09:02:18 2018 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Wed May 26 13:43:17 2021 (+0200)
##           By: Joerg Fallmann
##     Update #: 460
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

# Logging
import datetime
import glob

# multiprocessing
import multiprocessing

# numpy
import shlex
from itertools import repeat

# others
from natsort import natsorted
import pyBigWig as pbw

from RIssmed.Tweaks.logger import (
    makelogdir,
    makelogfile,
    listener_process,
    listener_configurer,
    worker_configurer,
)
from RIssmed import _version

__version__ = _version.get_versions()["version"]


# load own modules
from RIssmed.Tweaks.FileProcessor import *
from RIssmed.Tweaks.RNAtweaks import *
from RIssmed.Tweaks.RNAtweaks import _pl_to_array
from RIssmed.Tweaks.NPtweaks import *

log = logging.getLogger(__name__)  # use module name
SCRIPTNAME = os.path.basename(__file__).replace(".py", "")


def starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(apply_args_and_kwargs, args_for_starmap)


def apply_args_and_kwargs(fn, args, kwargs):
    return fn(*args, **kwargs)


def scan_input(
    queue,
    configurer,
    level,
    pat,
    cutoff,
    border,
    ulim,
    temperature,
    procs,
    unconstraint,
    unp,
    pai,
    outdir,
    dir,
    genes,
    chromsizes,
    padding,
):

    logid = SCRIPTNAME + ".scan_input: "
    try:
        # set path for output
        if outdir:
            log.info(logid + "Printing to " + outdir)
            if not os.path.isabs(outdir):
                outdir = os.path.abspath(outdir)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
        else:
            outdir = os.path.abspath(os.getcwd())

        pattern = pat.split(sep=",")
        window = int(pattern[0])
        span = int(pattern[1])

        genecoords = parse_annotation_bed_withchrom(
            genes
        )  # get genomic coords to print to bed later, should always be just one set of coords per gene

        log.debug(logid + str(genecoords))

        # Create process pool with processes
        num_processes = procs or 1
        call_list = []

        header = read_chromsize(chromsizes)
        rawbigfw = rawbigre = unpbigfw = unpbigre = paibigfw = paibigre = None

        if unconstraint:
            rawbigfw = pbw.open(os.path.join(outdir, f"{unconstraint}.fw.bw"), "w")
            rawbigfw.addHeader(header)
            rawbigre = pbw.open(os.path.join(outdir, f"{unconstraint}.re.bw"), "w")
            rawbigre.addHeader(header)
        if unp:
            unpbigfw = pbw.open(os.path.join(outdir, f"{unp}.fw.bw"), "w")
            unpbigfw.addHeader(header)
            unpbigre = pbw.open(os.path.join(outdir, f"{unp}.re.bw"), "w")
            unpbigre.addHeader(header)
        if pai:
            paibigfw = pbw.open(os.path.join(outdir, f"{pai}.fw.bw"), "w")
            paibigfw.addHeader(header)
            paibigre = pbw.open(os.path.join(outdir, f"{pai}.re.bw"), "w")
            paibigre.addHeader(header)

        for goi in genecoords:
            log.info(logid + "Working on " + goi)
            gs, ge, gstrand = get_location(genecoords[goi][0])

            raw = getfiles(unconstraint, window, span, temperature, goi)
            unpaired = getfiles("diffnu", window, span, temperature, goi)
            paired = getfiles("diffnp", window, span, temperature, goi)

            filelist = [raw, unpaired, paired]
            bwlist = [rawbigfw, rawbigre, unpbigfw, unpbigre, paibigfw, paibigre]

            call_list.append(
                (
                    goi,
                    filelist,
                    bwlist,
                    gs,
                    ge,
                    gstrand,
                    ulim,
                    cutoff,
                    border,
                    outdir,
                    padding,
                ),
            )

        with multiprocessing.Pool(num_processes, maxtasksperchild=1) as pool:
            outlist = starmap_with_kwargs(
                pool,
                generate_bws,
                call_list,
                repeat(
                    {
                        "queue": queue,
                        "configurer": configurer,
                        "level": level,
                    },
                    len(call_list),
                ),
            )
            pool.close()
            pool.join()
        for entry in outlist:
            writebws(entry, outdir, bwlist)

        for bw in [rawbigfw, rawbigre, unpbigfw, unpbigre, paibigfw, paibigre]:
            bw.close()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def generate_bws(
    goi,
    filelist,
    bwlist,
    gs,
    ge,
    gstrand,
    ulim,
    cutoff,
    border,
    outdir,
    padding,
    queue=None,
    configurer=None,
    level=None,
):

    logid = SCRIPTNAME + ".judge_diff: "
    try:
        if queue and level:
            configurer(queue, level)
        raw, up, pa = filelist

        if raw:
            chrom, strand, cons, reg, f, window, span, temperature = map(
                str, str(os.path.basename(raw)).replace(goi + "_", "", 1).split(sep="_")
            )
        elif up:
            chrom, strand, cons, reg, f, window, span, temperature = map(
                str, str(os.path.basename(up)).replace(f"StruCons_{goi}" + "_", "", 1).split(sep="_")
            )
        elif pa:
            chrom, strand, cons, reg, f, window, span, temperature = map(
                str, str(os.path.basename(pa)).replace(f"StruCons_{goi}" + "_", "", 1).split(sep="_")
            )

        span = span.split(sep=".")[0]
        cs, ce = map(int, cons.split(sep="-"))
        ws, we = map(int, reg.split(sep="-"))

        cs = cs - ws  # fit to window and make 0-based
        ce = ce - ws  # fit to window and make 0-based closed

        if 0 > any([cs, ce, ws, we]):
            raise Exception(
                "One of "
                + str([cs, ce, ws, we])
                + " lower than 0! this should not happen for "
                + ",".join([goi, chrom, strand, cons, reg, f, window, span, temperature])
            )

        if gstrand != "-":
            ws = ws + gs - 2  # get genomic coords 0 based closed, ws and gs are 1 based
            we = we + gs - 2

        else:
            wst = ws  # temp ws for we calc
            ws = ge - we  # get genomic coords 0 based closed, ge and we are 1 based
            we = ge - wst

        log.debug(
            logid
            + "Coords: "
            + " ".join(
                map(
                    str,
                    [
                        goi,
                        chrom,
                        strand,
                        cons,
                        reg,
                        f,
                        window,
                        span,
                        temperature,
                        gs,
                        ge,
                        cs,
                        ce,
                        ws,
                        we,
                    ],
                )
            )
        )

        border = abs(border)  # defines how big a diff has to be to be of importance

        log.info(
            logid
            + "Continuing "
            + str(goi)
            + " calculation with cutoff: "
            + str(cutoff)
            + " and border "
            + str(border)
        )  # + ' and ' + str(border2))

        # Read in plfold output
        noc = _pl_to_array(raw, ulim)

        out = {}
        out["r"] = {"fw": dict(), "re": dict()}
        out["u"] = {"fw": dict(), "re": dict()}
        out["p"] = {"fw": dict(), "re": dict()}

        if abs(np.nanmean(noc[cs : ce + 1])) <= cutoff:
            if up:
                uc = _pl_to_array(up, ulim)  # This is the diffacc for unpaired constraint
                # Calculate raw prop unpaired for constraint diff file
                uc = noc + uc
            if pa:
                pc = _pl_to_array(pa, ulim)  # This is the diffacc for paired constraint
                # Calculate raw prop unpaired for constraint diff file
                pc = noc + pc

            log.debug(
                logid
                + "unpaired: "
                + str(up)
                + " and paired: "
                + str(pa)
                + " Content: "
                + str(uc[ulim : ulim + 10])
                + " test "
                + str(np.all(uc[ulim : ulim + 10]))
            )

            """
            Collect positions of interest with padding around constraint
            Constraints are influencing close by positions strongest so strong influence of binding there is expected
            """

            for pos in range(len(noc)):
                if pos not in range(cs - padding + 1 - ulim, ce + padding + 1 + ulim):
                    if strand != "-":
                        gpos = pos + ws - ulim + 1  # already 0-based
                        gend = gpos + ulim  # 0-based half-open
                        orient = "fw"
                    else:
                        gpos = we - pos  # already 0-based
                        gend = gpos + ulim  # 0-based half-open
                        orient = "re"

                    if border < abs(noc[pos]):
                        out["raw"][orient]["chrom"].append(str(chrom))
                        out["raw"][orient]["start"].append(str(gpos))
                        out["raw"][orient]["end"].append(str(gend))
                        out["raw"][orient]["value"].append(noc[pos])

                    if border < abs(uc[pos]):
                        out["uc"][orient]["chrom"].append(str(chrom))
                        out["uc"][orient]["start"].append(str(gpos))
                        out["uc"][orient]["end"].append(str(gend))
                        out["uc"][orient]["value"].append(uc[pos])

                    if border < abs(uc[pos]):
                        out["pc"][orient]["chrom"].append(str(chrom))
                        out["pc"][orient]["start"].append(str(gpos))
                        out["pc"][orient]["end"].append(str(gend))
                        out["pc"][orient]["value"].append(pc[pos])

            return out

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def writebws(out, outdir, bwlist):

    logid = SCRIPTNAME + ".writebws: "

    try:
        if len(out["raw"]["fw"]["chrom"]) > 0:
            bw = bwlist[0]
            bw.addEntries(
                out["raw"]["fw"]["chrom"],
                out["raw"]["fw"]["start"],
                ends=out["raw"]["fw"]["end"],
                values=out["raw"]["fw"]["value"],
            )
        if len(out["raw"]["re"]["chrom"]) > 0:
            bw = bwlist[1]
            bw.addEntries(
                out["raw"]["re"]["chrom"],
                out["raw"]["re"]["start"],
                ends=out["raw"]["re"]["end"],
                values=out["raw"]["re"]["value"],
            )

        if len(out["uc"]["fw"]["chrom"]) > 0:
            bw = bwlist[2]
            bw.addEntries(
                out["uc"]["fw"]["chrom"],
                out["uc"]["fw"]["start"],
                ends=out["uc"]["fw"]["end"],
                values=out["uc"]["fw"]["value"],
            )
        if len(out["uc"]["re"]["chrom"]) > 0:
            bw = bwlist[3]
            bw.addEntries(
                out["uc"]["re"]["chrom"],
                out["uc"]["re"]["start"],
                ends=out["uc"]["re"]["end"],
                values=out["uc"]["re"]["value"],
            )

        if len(out["pc"]["fw"]["chrom"]) > 0:
            bw = bwlist[4]
            bw.addEntries(
                out["pc"]["fw"]["chrom"],
                out["pc"]["fw"]["start"],
                ends=out["pc"]["fw"]["end"],
                values=out["pc"]["fw"]["value"],
            )
        if len(out["pc"]["re"]["chrom"]) > 0:
            bw = bwlist[0]
            bw.addEntries(
                out["pc"]["re"]["chrom"],
                out["pc"]["re"]["start"],
                ends=out["pc"]["re"]["end"],
                values=out["pc"]["re"]["value"],
            )

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def getfiles(name, window, span, temperature, goi):
    """Retrieve files following specified pattern

    Parameters
    ----------
    name : str
        Name of file to search for, e.g. 'raw', 'unpaired' or 'paired'
    window : int
        size of folding window used
    span : int
        length of basepair span used
    temperature : float
        temperature used for folding
    goi : str
        Gene of interest

    Returns
    -------
    List
        List of files found folowing search pattern
    """

    ret = list()
    if not name:
        return None
    else:
        # get files with specified pattern
        lookfor = os.path.abspath(
            os.path.join(
                dir,
                goi,
                f"*{goi}*_{name}_*{str(window)}_{str(span)}_{str(temperature)}.npy",
            )
        )
        log.debug(logid + "PATHS: " + str(raw) + "\t" + str(paired) + "\t" + str(unpaired))
        collectall = natsorted(glob.glob(lookfor), key=lambda y: y.lower())
        # get absolute path for files
        fullname = [os.path.abspath(i) for i in collectall]
        log.debug(logid + "PATHS: " + str(len(r)) + "\t" + str(len(p)) + "\t" + str(len(u)))
        if not fullname:
            log.warning(
                logid
                + "Could not find files for Gene "
                + str(goi)
                + " and window "
                + str(window)
                + " and span "
                + str(span)
                + " and temperature "
                + str(temperature)
                + " Will skip"
            )
            return None
        else:
            return fullname


def read_chromsize(cs):
    logid = SCRIPTNAME + ".read_chromsize: "
    sizes = list()
    if ".gzip" in os.path.basename(cs)[-6:]:
        sizes = gzip.open(cs, "rt").read().splitlines()
    else:
        sizes = open(cs, "r").read().splitlines()
    sizes = [tuple((str(x), int(y))) for x, y in [l.split("\t") for l in sizes]]
    log.debug(f"{logid} {sizes}")
    return sizes


def main(args=None):
    """Main process, prepares run_settings dict, creates logging process queue and worker processes for folding, calls screen_genes

    Parameters
    ----------

    Returns
    -------
    Call to screen_genes
    """

    logid = SCRIPTNAME + ".main: "

    try:
        if not args:
            args = parseargs_browser()

        if args.version:
            sys.exit("Running RIssmed version " + __version__)

        #  Logging configuration
        logdir = args.logdir
        ts = str(datetime.datetime.now().strftime("%Y%m%d_%H_%M_%S_%f"))
        logfile = str.join(os.sep, [os.path.abspath(logdir), SCRIPTNAME + "_" + ts + ".log"])
        loglevel = args.loglevel

        makelogdir(logdir)
        makelogfile(logfile)

        queue = multiprocessing.Manager().Queue(-1)
        listener = multiprocessing.Process(
            target=listener_process,
            args=(queue, listener_configurer, logfile, loglevel),
        )
        listener.start()

        worker_configurer(queue, loglevel)

        log.info(logid + "Running " + SCRIPTNAME + " on " + str(args.procs) + " cores.")
        log.info(logid + "CLI: " + sys.argv[0] + " " + "{}".format(" ".join([shlex.quote(s) for s in sys.argv[1:]])))

        scan_input(
            queue,
            worker_configurer,
            loglevel,
            args.pattern,
            args.cutoff,
            args.border,
            args.ulimit,
            args.temperature,
            args.procs,
            args.unconstraint,
            args.unpaired,
            args.paired,
            args.outdir,
            args.dir,
            args.genes,
            args.chromsizes,
            args.padding,
        )

        queue.put_nowait(None)
        listener.join()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


####################
####    MAIN    ####
####################
if __name__ == "__main__":

    logid = SCRIPTNAME + ".main: "
    try:
        main()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))

######################################################################
### GenerateBigWig.py ends here
