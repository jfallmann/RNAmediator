#!/usr/bin/env python3
## CollectWindowResults.py ---
##
## Filename: CollectWindowResults.py
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
from itertools import repeat

# others
from natsort import natsorted

from RNAmediator.Tweaks.logger import (
    makelogdir,
    makelogfile,
    listener_process,
    listener_configurer,
    worker_configurer,
)

from RNAmediator import _version

__version__ = _version.get_versions()["version"]

# load own modules
from RNAmediator.Tweaks.FileProcessor import *
from RNAmediator.Tweaks.RNAtweaks import *
from RNAmediator.Tweaks.RNAtweaks import _get_ddg, _calc_ddg, _calc_gibbs, _calc_nrg
from RNAmediator.Tweaks.NPtweaks import *

log = logging.getLogger(__name__)  # use module name
SCRIPTNAME = os.path.basename(__file__).replace(".py", "")


def screen_genes(pat, border, procs, outdir, dir, genes, temperature, queue=None,
    configurer=None,
    level=None,):

    logid = SCRIPTNAME + ".screen_genes: "
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

        genecoords = parse_annotation_bed(
            genes
        )  # get genomic coords to print to bed later, should always be just one set of coords per gene

        # Create process pool with processes
        num_processes = procs or 1
        pool = multiprocessing.Pool(processes=num_processes, maxtasksperchild=1)

        for goi in genecoords:

            log.info(logid + "Working on " + goi)
            gs, ge, gstrand, value = get_location(genecoords[goi][0])

            # get files with specified pattern
            rawfiles = os.path.abspath(
                os.path.join(
                    dir,
                    goi,
                    goi + f"*{str(window)}_{str(span)}*_{str(temperature)}.gz",
                )
            )
            log.debug(f'COLLECTING: {rawfiles}')
            # search for files
            r = natsorted(glob.glob(rawfiles), key=lambda y: y.lower())
            log.debug(logid + "Files found: " + str(r))

            # get absolute path for files
            nocons = []

            rawfiles = [os.path.abspath(i) for i in r]

            if not rawfiles:
                log.warning(logid + "No output for gene " + str(goi) + " found, will skip!")
                continue

            try:
                for i in range(len(r)):
                    log.debug(logid + "Adding file " + str(r[i]) + " to queue.")
                    pool.apply_async(
                        calc,
                        args=(r[i], gs, ge, border, outdir),
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
       raise


def calc(p, gs, ge, border, outdir, queue=None, configurer=None, level=None):

    logid = SCRIPTNAME + ".calc_ddg: "
    try:
        log.debug(logid + "Parsing outputfile " + p)
        if queue and level:
            configurer(queue, level)

        goi, chrom, strand, cons, reg, window, span, temp = map(str, os.path.basename(p).split(sep="_"))
        border1, border2 = map(float, border.split(","))  # defines how big a diff has to be to be of importance

        log.info(logid + "Continuing calculation with borders: " + str(border1) + " and " + str(border2))

        out = defaultdict()
        ddgs = _get_ddg(p)
        log.debug(logid + "ddgs: " + str(ddgs))

        RT = (-1.9872041 * 10 ** (-3)) * (37 + 273.15)
        log.debug(logid + "RT is " + str(RT))

        cons = ddgs.pop("constraint")
        if not cons in out:
            out[cons] = list()
        ddg = _calc_ddg(ddgs)
        if ddg is not None:
            if ddg > border1 and ddg < border2:
                dkd = math.exp(ddg / RT)
                out[cons].append(
                    "\t".join(
                        [
                            str(chrom),
                            str(gs),
                            str(ge),
                            str(goi),
                            str(ddg),
                            str(strand),
                            str(cons),
                            str(dkd) + "\n",
                        ]
                    )
                )
        if out:
            write_out(out, outdir)
        else:
            log.warning(logid + "No ddg above cutoffs for gene " + str(goi))
        return

    except Exception:
       raise


def write_out(out, outdir):

    logid = SCRIPTNAME + ".savelist: "
    try:
        for cons in out:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            if not os.path.exists(os.path.abspath(os.path.join(outdir, "Collection_window.bed.gz"))):
                with gzip.open(
                    os.path.abspath(os.path.join(outdir, "Collection_window.bed.gz")),
                    "wb",
                ) as o:
                    o.write(bytes("\n".join(out[cons]), encoding="UTF-8"))
            else:
                with gzip.open(
                    os.path.abspath(os.path.join(outdir, "Collection_window.bed.gz")),
                    "ab",
                ) as o:
                    o.write(bytes("\n".join(out[cons]), encoding="UTF-8"))
    except Exception:
        raise


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
            args = parseargs_collect_window()

        if args.version:
            sys.exit("Running RNAmediator version " + __version__)

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

        try:
            result = screen_genes(
                args.pattern,
                args.border,
                args.procs,
                args.outdir,
                args.dir,
                args.genes,
                args.temperature,
                queue=queue,
                configurer=worker_configurer,
                level=args.loglevel,
            )
            queue.put(None)
            listener.join()
            return 0
        except Exception:
            raise

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        print(f'ERROR: {logid} {"".join(tbe.format())}')
        if log:
            log.error(logid + "".join(tbe.format()))
        listener.terminate()
        sys.exit(1)
        

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
