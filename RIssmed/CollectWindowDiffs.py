#!/usr/bin/env python3
## CollectWindowDiffs.py ---
##
## Filename: CollectWindowDiffs.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Aug 15 13:49:46 2019 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Mon May 11 15:46:56 2020 (+0200)
##           By: Joerg Fallmann
##     Update #: 72
## URL:
## Doc URL:
## Keywords:
## Compatibility:
##
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
import itertools

# others
from natsort import natsorted

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
from RIssmed.Tweaks.RNAtweaks import _get_ddg, _calc_ddg, _calc_gibbs, _calc_nrg, _pl_to_array
from RIssmed.Tweaks.NPtweaks import *

log = logging.getLogger(__name__)  # use module name
SCRIPTNAME = os.path.basename(__file__).replace(".py", "")


def screen_diffs(queue, configurer, level, window, span, ulim, cutoff, procs, outdir, genes):

    logid = SCRIPTNAME + ".screen_diffs: "
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

        windows = window.split(",")
        spans = span.split(",")

        combis = [x for x in [itertools.product(w + "_", spans) for w in windows]]
        patterns = itertools.combinations(combis, 2)

        genecoords = parse_annotation_bed(
            genes
        )  # get genomic coords to print to bed later, should always be just one set of coords per gene

        # Create process pool with processes
        num_processes = procs or 1
        pool = multiprocessing.Pool(processes=num_processes, maxtasksperchild=1)

        for goi in genecoords:
            log.info(logid + "Working on " + goi)
            gs, ge, gstrand = get_location(genecoords[goi][0])

            # get files with specified pattern
            for pair in patterns:
                diffin = [
                    os.path.abspath(
                        os.path.join(
                            goi,
                            goi + "*" + x + ".npy",
                        )
                    )
                    for x in pair
                ]

                # search for files
                p = [natsorted(glob.glob(d), key=lambda y: y.lower()) for d in diffin]
                log.debug(logid + "Files found: " + str(p))

                # get absolute path for files

                comparelist = [os.path.abspath(i) for i in p]

                if not comparelist:
                    log.warning(logid + "No output for gene " + str(goi) + " found, will skip!")
                    continue

                try:
                    pool.apply_async(
                        calcdiff,
                        args=(comparelist, gs, ge, gstrand, ulim, cutoff, outdir),
                        kwds={"queue": queue, "configurer": configurer, "level": level},
                    )
                except Exception:
                    exc_type, exc_value, exc_tb = sys.exc_info()
                    tbe = tb.TracebackException(
                        exc_type,
                        exc_value,
                        exc_tb,
                    )
                    log.error(logid + "".join(tbe.format()))

        pool.close()
        pool.join()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def calcdiff(p, gs, ge, gstrand, ulim, border, outdir, queue=None, configurer=None, level=None):

    logid = SCRIPTNAME + ".calc_diff: "
    try:
        log.debug(logid + "Parsing outputfile " + p)
        if queue and level:
            configurer(queue, level)

        goi, chrom, strand, cons, reg, f, window, span = map(str, os.path.basename(p[0]).split(sep="_"))
        goi1, chrom1, strand1, cons1, reg1, f1, window1, span1 = map(str, os.path.basename(p[1]).split(sep="_"))
        span = span.split(sep=".")[0]
        span1 = span1.split(sep=".")[0]
        diffname = "_".join(["|".join([window, window1]), "|".join([span, span1])])

        cs, ce = map(int, cons.split(sep="-"))
        ws, we = map(int, reg.split(sep="-"))

        cs = cs - ws  # fit to window and make 0-based
        ce = ce - ws  # fit to window and make 0-based closed

        if 0 > any([cs, ce, ws, we]):
            raise Exception(
                "One of "
                + str([cs, ce, ws, we])
                + " lower than 0! this should not happen for "
                + ",".join([goi, chrom, strand, cons, reg, f, window, span])
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
            + "DiffCoords: "
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
            logid + "Continuing " + str(goi) + " calculation with cutoff: " + str(border)
        )  # + ' and ' + str(border2))

        out = defaultdict()
        a1 = _pl_to_array(p[0], ulim)
        a2 = _pl_to_array(p[1], ulim)
        if not np.array_equal(a1, a2):
            diff = a1 - a2
        else:
            log.info(logid + "No influence on structure with unpaired constraint at " + cons)
            diff = None

        log.debug(logid + "diff: " + str(diff[1:10]))

        # RT = (-1.9872041 * 10 ** (-3)) * (37 + 273.15)
        # log.debug(logid + "RT is " + str(RT))

        if diff is not None:
            for pos in range(len(diff)):
                if strand != "-":
                    gpos = pos + ws - ulim + 1  # already 0-based
                    gend = gpos + ulim  # 0-based half-open
                    gcst = cs + ws + 1
                    gcen = ce + ws + 2
                    gcons = str(gcst) + "-" + str(gcen)
                else:
                    gpos = we - pos  # already 0-based
                    gend = gpos + ulim  # 0-based half-open
                    gcst = we - ce - 1
                    gcen = we - cs
                    gcons = str(gcst) + "-" + str(gcen)

                if border < abs(diff[pos]):
                    out[diffname].append(
                        "\t".join(
                            [
                                str(chrom),
                                str(gpos),
                                str(gend),
                                str(goi) + "|" + str(cons) + "|" + str(gcons),
                                str(diff[pos]),
                                str(strand) + "\n",
                            ]
                        )
                    )
        if len(out) > 0:
            write_out(out, outdir)
        else:
            log.warning(logid + "No diffs above cutoffs for gene " + str(goi))
        return
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def write_out(out, outdir):

    logid = SCRIPTNAME + ".savelist: "
    try:
        for cons in out:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            if not os.path.exists(os.path.abspath(os.path.join(outdir, "Collection_windowdiff.bed.gz"))):
                with gzip.open(
                    os.path.abspath(os.path.join(outdir, "Collection_windowdiff.bed.gz")),
                    "wb",
                ) as o:
                    o.write(bytes("\n".join(out[cons]), encoding="UTF-8"))
            else:
                with gzip.open(
                    os.path.abspath(os.path.join(outdir, "Collection_windowdiff.bed.gz")),
                    "ab",
                ) as o:
                    o.write(bytes("\n".join(out[cons]), encoding="UTF-8"))
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


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
            args = parseargs_collect_windowdiff()

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

        screen_diffs(
            queue,
            worker_configurer,
            loglevel,
            args.window,
            args.span,
            args.ulimit,
            args.cutoff,
            args.procs,
            args.outdir,
            args.genes,
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
#
# CollectWindowResults.py ends here
