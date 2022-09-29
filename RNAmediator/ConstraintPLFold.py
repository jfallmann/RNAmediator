#!/usr/bin/env python3
### ConstraintPLFold.py ---
##
## Filename: ConstraintPLFold.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Sep  6 09:02:18 2018 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Sat Jun 25 15:14:04 2022 (+0200)
##           By: Joerg Fallmann
##     Update #: 449
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
############
# find ffmpeg executable
# import shutil
# plt.rcParams['animation.ffmpeg_path'] = shutil.which("ffmpeg")
# plt.rc('verbose', level='debug-annoying', fileo=sys.stdout)
# matplotTweaks.verbose.set_level("helpful")
# plt.rc('animation', html='html5')
##
######################################################################
##
### Change Log:
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
## IMPORTS
from __future__ import annotations  # It will become the default in Python 3.10

# numpy
# RNA
# Logging
import errno
import multiprocessing
import shlex
import re
from typing import Dict
from RNAmediator import _version

__version__ = _version.get_versions()["version"]

# load own modules
from RNAmediator.Tweaks.FileProcessor import *
from RNAmediator.Tweaks.RNAmediator import (
    preprocess,
    SequenceSettings,
    rnamediator_logging_setup,
    expand_pl_window,
    localize_pl_window,
)
from RNAmediator.Tweaks.RNAtweaks import *
from RNAmediator.Tweaks.RNAtweaks import _npprint
from RNAmediator.Tweaks.NPtweaks import *

# Biopython stuff

log = logging.getLogger(__name__)  # use module name
# log.propagate = True
SCRIPTNAME = os.path.basename(__file__).replace(".py", "")


def pl_fold(
    window,
    span,
    region,
    temperature,
    multi,
    unconstrained,
    unpaired,
    paired,
    save,
    procs,
    outdir,
    run_settings: Dict[str, SequenceSettings],
    pattern=None,
    queue=None,
    configurer=None,
    level=None,
):
    """Run RNAplfold

    Parameters
    ----------
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    region: int -u Option of RNAplfold
    temp: float Temperature
    multi: int Multiplyer for Window extension
    unconstrained: str Infix for ouput files of unconstrained folding
    unpaired: str Infix for ouput files of unpaired constraint
    paired: str Infix for output files of paired constraint, None will skip paired constraint folding
    save: bool Should npy or npy and gz output be saved
    procs: int Number of cores to use
    outdir : str Location of the Output directory. If it is an empty string os.cwd() is used
    run_settings: Dict[str, SequenceSettings] Dictionary with sequence settings per Fasta ID
    pattern: str Filter for Gene of Interest (optional)
    queue: multiprocessing_queue Logging process queue
    configurer: multiprocessing_config for Logging processes
    level: logging.level Level for log process

    """
    logid = SCRIPTNAME + ".plfold: "
    try:
        log.debug(f"q:{queue} c:{configurer} l:{level} r:{run_settings}")
        if queue and level:
            configurer(queue, level)
        # Create process pool with processes
        num_processes = procs or 1
        # with get_context("spawn").Pool(processes=num_processes-1, maxtasksperchild=1) as pool:
        pool = multiprocessing.Pool(processes=num_processes, maxtasksperchild=1)

        # Start the work
        for fasta in run_settings:
            fasta_settings = run_settings[fasta]
            goi = fasta_settings.gene
            gs = fasta_settings.genomic_coords.start
            ge = fasta_settings.genomic_coords.end
            gstrand = fasta_settings.genomic_coords.strand
            seq_record = fasta_settings.sequence_record
            conslist = fasta_settings.constrainlist
            log.debug(logid + str(conslist))

            if len(seq_record.seq) < window * multi:
                log.warning(
                    str(
                        "Sequence of "
                        + goi
                        + " too short, seqlenght "
                        + str(len(seq_record.seq))
                        + " with window size "
                        + str(window)
                        + " and multiplyer "
                        + str(multi)
                    )
                )
                # continue

            if pattern and pattern not in goi:
                continue

            log.info(f"{logid} Working on {goi}: {seq_record.id}, start {gs}, end {ge}, strand {gstrand} ")

            # define data structures
            data = {"up": []}
            an = [np.nan]
            # We check if we need to fold the whole seq or just a region around the constraints
            conslist = fasta_settings.constrainlist
            constype = fasta_settings.constraintype
            log.debug(f"{logid} {conslist} {constype}")

            for cons_tuple in conslist:
                log.debug(logid + "ENTRY: " + str(cons_tuple))
                if cons_tuple[0] is None:  # in case we just want to fold the sequence without constraints at all
                    # raise NotImplementedError("Needs to be reimplemented")
                    cons = str(cons_tuple[0])
                    log.info(f"{logid} Folding without constraint from {gs} to {ge}")
                    start, end = [1, ge - gs]
                    tostart, toend = expand_pl_window(start, end, window, multi, len(seq_record.seq))
                    cons = str(start) + "-" + str(end) + "_" + str(tostart) + "-" + str(toend)
                    log.debug(logid + str.join(" ", [goi, cons, gstrand]))

                    if start < 1 or end > len(seq_record.seq):
                        log.warning(
                            logid
                            + "Constraint out of sequence bounds! skipping! "
                            + ",".join(
                                map(
                                    str,
                                    [
                                        goi,
                                        len(seq_record.seq),
                                        str(start) + "-" + str(end),
                                    ],
                                )
                            )
                        )
                        continue
                    if check_raw_existing(
                        str(seq_record.id),
                        unconstrained,
                        cons,
                        region,
                        window,
                        span,
                        temperature,
                        outdir,
                    ):
                        log.warning(logid + str(cons) + " Exists for " + str(seq_record.id) + "! Skipping!")
                        continue
                    pool.apply_async(
                        scan_seq,
                        args=(
                            seq_record.id,
                            seq_record.seq,
                            start,
                            end,
                            window,
                            span,
                            region,
                            temperature,
                            multi,
                            save,
                            outdir,
                        ),
                        kwds={
                            "unconstrained": unconstrained,
                            "queue": queue,
                            "configurer": configurer,
                            "level": level,
                        },
                    )

                else:
                    if len(cons_tuple) > 1:  # Then paired constraints should be used
                        cons_tuple = [str(cons) for cons in cons_tuple]
                        cons = ":".join(cons_tuple)
                        fstart, fend = [None, None]
                        start, end = [None, None]

                        if gstrand == "+" or gstrand == ".":
                            [fstart, fend], [start, end] = [
                                [x - gs for x in get_location(cn)[:2]] for cn in cons.split(":", 1)
                            ]
                        else:
                            [fstart, fend], [start, end] = [
                                [ge - x for x in get_location(cn)[:2][::-1]] for cn in cons.split(":", 1)
                            ]
                        consval = str([get_location(cn)[3] for cn in cons.split(":", 1)])
                        cons = str(fstart) + "-" + str(fend) + ":" + str(start) + "-" + str(end)
                        if start < 0 or fstart < 0 or end > len(seq_record.seq) or fend > len(seq_record.seq):
                            log.warning(
                                logid
                                + "Constraint out of sequence bounds! skipping! "
                                + ",".join(
                                    map(
                                        str,
                                        [
                                            goi,
                                            len(seq_record.seq),
                                            str(start) + "-" + str(end),
                                            str(fstart) + "-" + str(fend),
                                        ],
                                    )
                                )
                            )
                            continue
                        if checkexisting(
                            str(seq_record.id),
                            paired,
                            unpaired,
                            cons,
                            region,
                            window,
                            span,
                            temperature,
                            outdir,
                        ):
                            log.warning(logid + str(cons) + " Exists for " + str(seq_record.id) + "! Skipping!")
                            continue

                        log.info(logid + "Constraining to " + str(fstart) + " and " + str(fend))
                        goi, chrom, strand = idfromfa(seq_record.id)
                        pool.apply_async(
                            constrain_seq_paired,
                            args=(
                                seq_record.id,
                                seq_record.seq,
                                fstart,
                                fend,
                                start,
                                end,
                                window,
                                span,
                                region,
                                temperature,
                                multi,
                                paired,
                                unpaired,
                                save,
                                outdir,
                                constype,
                                consval,
                            ),
                            kwds={
                                "unconstrained": unconstrained,
                                "queue": queue,
                                "configurer": configurer,
                                "level": level,
                            },
                        )

                    else:
                        # indexing because conslist is a list of tuples for multi constraints
                        cons = str(cons_tuple[0])
                        log.info(logid + "Calculating constraint\t" + cons)
                        if gstrand == "+" or gstrand == ".":
                            start, end = [x - gs for x in get_location(cons)[:2]]
                        else:
                            start, end = [ge - x for x in get_location(cons)[:2][::-1]]

                        tostart, toend = expand_pl_window(start, end, window, multi, len(seq_record.seq))
                        consval = str([get_location(cn)[3] for cn in cons.split(":", 1)])
                        cons = str(start) + "-" + str(end) + "_" + str(tostart) + "-" + str(toend)
                        log.debug(logid + str.join(" ", [goi, cons, gstrand]))

                        if start < 0 or end > len(seq_record.seq):
                            log.warning(
                                logid
                                + "Constraint out of sequence bounds! skipping! "
                                + ",".join(
                                    map(
                                        str,
                                        [
                                            goi,
                                            len(seq_record.seq),
                                            str(start) + "-" + str(end),
                                        ],
                                    )
                                )
                            )
                            continue
                        if checkexisting(
                            str(seq_record.id),
                            paired,
                            unpaired,
                            cons,
                            region,
                            window,
                            span,
                            temperature,
                            outdir,
                        ):
                            log.warning(logid + str(cons) + " Exists for " + str(seq_record.id) + "! Skipping!")
                            continue

                        pool.apply_async(
                            constrain_seq,
                            args=(
                                seq_record.id,
                                seq_record.seq,
                                start,
                                end,
                                window,
                                span,
                                region,
                                temperature,
                                multi,
                                paired,
                                unpaired,
                                save,
                                outdir,
                                constype,
                                consval,
                            ),
                            kwds={
                                "unconstrained": unconstrained,
                                "queue": queue,
                                "configurer": configurer,
                                "level": level,
                            },
                        )

        pool.close()        
        pool.join()  # timeout
        log.info(logid + "DONE: output in: " + str(outdir))

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))
        pool.terminate()
        pool.join()
        queue.terminate()
        queue.join()
        sys.exit(1)

    return 0


def fold_unconstraint(
    seq,
    id,
    region,
    window,
    span,
    temperature,
    unconstrained,
    save,
    outdir,
    rawentry=None,
    locws=None,
    locwe=None,
    queue=None,
    configurer=None,
    level=None,
):
    """Fold Sequence without constraint

    Parameters
    ----------
    seq: str Sequence to fold
    id: str Id for Sequence
    region: int -u Option of RNAplfold
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    temperature: float Temperature
    unconstrained: str Infix for ouput files of unconstrained folding
    save: bool Should npy or npy and gz output be saved
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    rawentry: str Name for raw file to pass to write_unconstraint
    locws: int Local Window start
    locwe: int Local Window end
    queue: multiprocessing_queue Logging process queue
    configurer: multiprocessing_config for Logging processes
    level: logging.level Level for log process

    Returns
    -------
    PLFoldOutput: PLFoldOutput object
    """
    logid = SCRIPTNAME + ".fold_unconstraint: "
    try:
        seq = seq.upper().replace("T", "U")
        if queue and level:
            configurer(queue, level)

        log.debug(
            f"{SCRIPTNAME} Seqlength: {len(seq)}, Window: {window}, Span: {span}, Region: {region}, Temperature: {temperature}"
        )
        plfold_output = api_rnaplfold(seq, window, span, region, temperature)

        if locws is not None and locwe is not None:  # If we only need a subset of the folded sequence
            log.debug(logid + "Cutting RIO from fold with boundaries " + str(locws) + " and " + str(locwe))
            plfold_output.localize(locws, locwe + 1)
            seq = seq[locws - 1 : locwe]

        write_unconstraint(
            save,
            str(id),
            str(seq),
            unconstrained,
            plfold_output,
            int(region),
            str(window),
            str(span),
            str(temperature),
            outdir,
            rawentry,
        )
        return plfold_output

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def scan_seq(
    sid,
    seq,
    start,
    end,
    window,
    span,
    region,
    temperature,
    multi,
    save,
    outdir,
    unconstrained=None,
    queue=None,
    configurer=None,
    level=None,
):
    """Fold Sequences without constraint for structure profile generation

    Parameters
    ----------
    sid: str Id for Sequence
    seq: str Sequence to fold
    start: int Start of subsequence to fold
    end: int End of subsequence to fold
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    region: int -u Option of RNAplfold
    temperature: float Temperature
    multi: int Multiplyer for window extension
    save: bool Should npy or npy and gz output be saved
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    unconstrained: str Infix for ouput files of unconstrained folding
    queue: multiprocessing_queue Logging process queue
    configurer: multiprocessing_config for Logging processes
    level: logging.level Level for log process

    Returns
    -------
    Call to fold_unconstraint, 1 on success
    """
    sid = str(sid)
    seq = str(seq)
    logid = SCRIPTNAME + ".scan_seq: "

    try:
        if queue and level:
            configurer(queue, level)

        seq = seq.upper().replace("T", "U")
        goi, chrom, strand = idfromfa(sid)
        log.debug(logid + "Scanning Sequence")

        # Here we fold the whole raw sequence +- window size

        tostart, toend = expand_pl_window(start, end, window, multi, len(seq))
        log.debug(f"{logid} tostart {tostart}, toend {toend}")
        seqtofold = str(seq[tostart - 1 : toend])
        log.debug(f"{logid}: Seq: {len(seq)}, fold: {len(seqtofold)}")
        # get local window of interest 0 based closed, we do not need to store the whole seqtofold
        locws, locwe = localize_pl_window(start, end, window, len(seq))
        cons = str("-".join([str(start), str(end)]) + "_" + "-".join([str(locws), str(locwe)]))

        if len(seqtofold) < (toend - tostart):
            log.warning(
                logid + "Sequence to small, skipping " + str(sid) + "\t" + str(len(seqtofold)) + "\t" + str(cons)
            )
            return

        log.debug(logid + str.join(" ", [goi, cons, strand]))

        if start < 1 or end > len(seq):
            log.warning(logid + f"Constraint out of sequence bounds! skipping! {len(seq)}, {str(start)} - {str(end)}")
            return

        if check_raw_existing(sid, unconstrained, cons, region, window, span, temperature, outdir):
            log.warning(logid + str(cons) + " Existst for " + str(sid) + "! Skipping!")
            return

        # get local start, ends 0 based closed
        locstart = start - tostart
        locend = end - tostart

        log.debug(
            " ".join(
                map(
                    str,
                    [
                        logid,
                        sid,
                        region,
                        str(len(seq)),
                        str(len(seqtofold)),
                        cons,
                        tostart,
                        locstart,
                        locend,
                        toend,
                    ],
                )
            )
        )

        # Cut sequence of interest from data, we no longer need the window extension as no effect outside of window
        # is visible with plfold anyways
        # get local start,ends 0 based closed
        locws = locws - tostart
        locwe = locwe - tostart

        plfold_unconstraint = fold_unconstraint(
            str(seq), sid, region, window, span, temperature, unconstrained, save, outdir
        )

        return 1

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def constrain_seq(
    sid,
    seq,
    start,
    end,
    window,
    span,
    region,
    temperature,
    multi,
    paired,
    unpaired,
    save,
    outdir,
    constype,
    consval,
    unconstrained=None,
    queue=None,
    configurer=None,
    level=None,
):
    """Fold Sequences with constraint

    Parameters
    ----------
    sid: str Id for Sequence
    seq: str Sequence to fold
    start: int Start of subsequence to fold
    end: int End of subsequence to fold
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    region: int -u Option of RNAplfold
    temperature: float Temperature
    multi: int Multiplyer for window extension
    paired: str Infix for output files of paired constraint
    unpaired: str Infix for ouput files of unpaired constraint
    save: bool Should npy or npy and gz output be saved
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    constype : str
        Type of constraint to apply, can be ['hard' (Default), 'soft'(not implemented yet), 'mutate']
    consval : str
        Value for constraint to apply, ignored when constype is 'hard'
    unconstrained: str, optional
        Infix for ouput files of unconstrained folding
    queue: multiprocessing_queue, optional
        Logging process queue
    configurer: multiprocessing_config, optional
        Config for Logging processes
    level: logging.level, optional
        Level for log process

    Returns
    -------
    Call to write_constraint
    """

    logid = SCRIPTNAME + ".constrain_seq: "
    seq = str(seq)
    sid = str(sid)
    try:
        if queue and level:
            configurer(queue, level)

        seq = seq.upper().replace("T", "U")
        goi, chrom, strand = idfromfa(sid)
        log.debug(logid + "CONSTRAINING AWAY with " + str(start) + " " + str(end))

        # for all constraints we now extract subsequences to compare against
        # we no longer fold the whole raw sequence but only the constraint region +- window size
        # Yes this is a duplicate but not if used in other context as standalone function
        tostart, toend = expand_pl_window(start, end, window, multi, len(seq))
        seqtofold = str(seq[tostart - 1 : toend])

        # get local window of interest 0 based closed, we do not need to store the whole seqtofold
        locws, locwe = localize_pl_window(start, end, window, len(seq))
        cons = str("-".join([str(start), str(end)]) + "_" + "-".join([str(locws), str(locwe)]))

        if len(seqtofold) < (toend - tostart):
            log.warning(
                logid + "Sequence to small, skipping " + str(sid) + "\t" + str(len(seqtofold)) + "\t" + str(cons)
            )
            return

        log.debug(logid + str.join(" ", [goi, cons, strand]))

        if start < 1 or end > len(seq):
            log.warning(
                logid
                + "Constraint out of sequence bounds! skipping! "
                + ",".join([str(len(seq)), str(start) + "-" + str(end)])
            )
            return

        if checkexisting(sid, paired, unpaired, cons, region, window, span, temperature, outdir):
            log.warning(logid + str(cons) + " Existst for " + str(sid) + "! Skipping!")
            return

        # get local start,ends 0 based closed
        locstart = start - tostart
        locend = end - tostart
        value = consval if constype != "hard" else None

        log.debug(
            " ".join(
                map(
                    str,
                    [
                        logid,
                        sid,
                        region,
                        str(len(seq)),
                        str(len(seqtofold)),
                        cons,
                        tostart,
                        locstart,
                        locend,
                        toend,
                        value,
                        unpaired,
                        paired,
                    ],
                )
            )
        )

        # Cut sequence of interest from data, we no longer need the window extension as no effect outside of window
        # is visible with plfold anyways
        # get local start,ends 0 based closed
        locws = locws - tostart
        locwe = locwe - tostart

        plfold_unconstraint = fold_unconstraint(
            str(seqtofold),
            sid,
            region,
            window,
            span,
            temperature,
            unconstrained,
            save,
            outdir,
            cons,
            locws,
            locwe,
        )
        an = plfold_unconstraint.get_rnamediator_np_array()  # create numpy array from output

        # fold constrained
        plfold_unpaired = api_rnaplfold(
            seqtofold,
            window,
            span,
            region,
            temperature,
            constype,
            consval,
            constraint=[("unpaired", locstart, locend + 1)],
        )
        plfold_unpaired.localize(locws, locwe + 1)
        au = plfold_unpaired.get_rnamediator_np_array()

        if paired:
            if constype in ["hard", "mutate"]:
                plfold_paired = api_rnaplfold(
                    seqtofold,
                    window,
                    span,
                    region,
                    temperature,
                    constype,
                    consval,
                    constraint=[("paired", locstart, locend + 1)],
                )
                plfold_paired.localize(locws, locwe + 1)
                ap = plfold_paired.get_rnamediator_np_array()
            else:
                raise NotImplementedError("Soft constraints for paired sequences need specific basepairs")
        else:
            plfold_paired = None
            ap = np.empty(au.shape)

        # Calculating accessibility difference between unconstrained and constraint fold, <0 means less accessible
        # with constraint, >0 means more accessible upon constraint

        diff_nu = diff_np = None

        if not np.array_equal(an, au):
            diff_nu = au - an
        else:
            log.info(logid + "No influence on Structure with unpaired constraint at " + cons)

        if constype in ["hard", "mutate"] and paired:
            if not np.array_equal(an, ap):
                diff_np = ap - an
            else:
                log.info(logid + "No influence on Structure with paired constraint at " + cons)

        seqtoprint = seqtofold[locws - 1 : locwe]

        write_constraint(
            save,
            str(sid),
            seqtoprint,
            paired,
            unpaired,
            plfold_unpaired,
            plfold_paired,
            cons,
            int(region),
            diff_nu,
            diff_np,
            str(window),
            str(span),
            str(temperature),
            outdir,
        )

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def constrain_seq_paired(
    sid,
    seq,
    fstart,
    fend,
    start,
    end,
    window,
    span,
    region,
    temperature,
    multi,
    paired,
    unpaired,
    save,
    outdir,
    constype,
    consval,
    unconstrained=None,
    queue=None,
    configurer=None,
    level=None,
):
    """Fold Sequences with paired constraint

    Parameters
    ----------
    sid: str Id for Sequence
    seq: str Sequence to fold
    fstart: int Start of First subsequence to fold
    fend: int End of First subsequence to fold
    start: int Start of subsequence to fold
    end: int End of subsequence to fold
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    region: int -u Option of RNAplfold
    temperature: float Temperature
    multi: int Multiplyer for window extension
    paired: str Infix for output files of paired constraint
    unpaired: str Infix for ouput files of unpaired constraint
    save: bool Should npy or npy and gz output be saved
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    constype : str
        Type of constraint to apply, can be ['hard' (Default), 'soft'(not implemented yet), 'mutate']
    consval : str
        Value for constraint to apply, ignored when constype is 'hard'
    unconstrained: str, optional
        Infix for ouput files of unconstrained folding
    queue: multiprocessing_queue, optional
        Logging process queue
    configurer: multiprocessing_config, optional
        Config for Logging processes
    level: logging.level, optional
        Level for log process
    Returns
    -------
    Call to write_constraint
    """

    logid = SCRIPTNAME + ".constrain_seq_paired: "
    seq = str(seq)
    sid = str(sid)
    try:
        if queue and level:
            configurer(queue, level)

        seq = str(seq).upper().replace("T", "U")
        log.debug(
            f"{logid} i:{sid} s:{seq} f:{fstart} e:{fend} s:{start} e:{end} w:{window} s:{span} u:{region} t:{temperature} m:{multi} p:{paired} u:{unpaired} a:{save} o:{outdir} c:{constype} v:{consval}"
        )
        goi, chrom, strand = idfromfa(sid)
        log.debug(logid + "CONSTRAINING PAIRWISE with " + str(start) + " " + str(end))

        # we no longer fold the whole sequence but only the constraint region +- window size
        tostart, toend = expand_pl_window(start, fend, window, multi, len(seq))
        seqtofold = str(seq[tostart - 1 : toend])

        # get local window of interest 0 based closed, we do not need to store the whole seqtofold
        locws, locwe = localize_pl_window(start, fend, window, len(seq))
        cons = str(
            "-".join([str(start), str(end) + ":" + str(fstart), str(fend)]) + "_" + "-".join([str(locws), str(locwe)])
        )

        if len(seqtofold) < (toend - tostart):
            log.warning(
                logid + "Sequence to small, skipping " + str(sid) + "\t" + str(len(seqtofold)) + "\t" + str(cons)
            )
            return

        if start < 1 or end > len(seq) or fstart < 1 or fend > len(seq):
            log.warning(
                logid
                + "Constraint out of sequence bounds! skipping! "
                + ",".join(
                    [
                        len(seq),
                        str(start) + "-" + str(end),
                        str(fstart) + "-" + str(fend),
                    ]
                )
            )
            return

        if checkexisting(sid, paired, unpaired, cons, region, window, span, temperature, outdir):
            log.warning(logid + str(cons) + " Exists for " + str(sid) + "! Skipping!")
            return

        # get local start,ends 0 based closed
        locstart = start - tostart
        locend = end - tostart
        flocstart = fstart - tostart
        flocend = fend - tostart
        value = consval if constype != "hard" else None

        log.debug(
            " ".join(
                map(
                    str,
                    [
                        logid,
                        sid,
                        region,
                        str(len(seq)),
                        str(len(seqtofold)),
                        cons,
                        tostart,
                        locstart,
                        locend,
                        toend,
                        value,
                    ],
                )
            )
        )

        locws = locws - tostart
        locwe = locwe - tostart

        # fold unconstrained
        plfold_unconstraint = fold_unconstraint(
            str(seqtofold),
            sid,
            region,
            window,
            span,
            temperature,
            unconstrained,
            save,
            outdir,
            cons,
            locws,
            locwe,
        )
        an = plfold_unconstraint.get_rnamediator_np_array()

        # fold constrained
        plfold_unpaired = api_rnaplfold(
            seqtofold,
            window,
            span,
            region,
            temperature,
            constype,
            consval,
            constraint=[
                ("unpaired", flocstart, flocend + 1),
                ("unpaired", locstart, locend + 1),
            ],
        )
        plfold_unpaired.localize(locws, locwe + 1)
        au = plfold_unpaired.get_rnamediator_np_array()

        if constype in ["hard", "mutate"] and paired:
            plfold_paired = api_rnaplfold(
                seqtofold,
                window,
                span,
                region,
                temperature,
                constype,
                consval,
                constraint=[
                    ("paired", flocstart, flocend + 1),
                    ("paired", locstart, locend + 1),
                ],
            )
            plfold_paired.localize(locws, locwe + 1)
            ap = plfold_paired.get_rnamediator_np_array()
        else:
            plfold_paired = None
        # Calculating accessibility difference between unconstrained and constraint fold, <0 means less accessible with constraint, >0 means more accessible upon constraint

        if not np.array_equal(an, au):
            diff_nu = au - an
        else:
            log.info(logid + "No influence on Structure with unpaired constraint at " + cons)
            diff_nu = None
        if not np.array_equal(an, ap) and constype in ["hard", "mutate"]:
            diff_np = ap - an
        else:
            if constype in ["hard", "mutate"]:
                log.info(logid + "No influence on Structure with paired constraint at " + cons)
            diff_np = None

        seqtoprint = seqtofold[locws - 1 : locwe]

        write_constraint(
            save,
            str(sid),
            seqtofold,
            paired,
            unpaired,
            plfold_unpaired,
            plfold_paired,
            cons,
            int(region),
            diff_nu,
            diff_np,
            str(window),
            str(span),
            str(temperature),
            outdir,
        )

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


# def constrain_temp(sid, seq, temp, window, span, region, multi, an, save, outdir, queue=None, configurer=None, level=None):
#     seq = seq.upper().replace("T", "U")
#     logid = SCRIPTNAME + '.constrain_temp: '
#     try:
#         if queue and level:
#             configurer(queue, level)
#
#         if len(seq) < int(window):
#             log.warning(logid+'Sequence to small, skipping '+str(sid)+'\t'+str(temp))
#             return
#
#         #refresh model details
#         #RNA = importTweaks.import_module('RNA')
#         md = RNA.md()
#         md.max_bp_span = span
#         md.window_size = window
#         #set temperature
#         md.temperature = temp
#         #create new fold_compound objects
#         fc_t = RNA.fold_compound(str(seq), md, RNA.OPTION_WINDOW)
#         #new data struct
#         data_t = {'up': []}
#
#         fc_t.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_t)
#
#         at = up_to_array(data_t['up'],region,len(seq))
#
#         diff_nt = an - at
#
#         #####Set temp, and rewrite save for temp
#         write_temp(save, sid, seq, str(temp), data_t, int(region), diff_nt, str(window), str(span), outdir)
#
#         return 1
#
#     except Exception:
#         exc_type, exc_value, exc_tb = sys.exc_info()
#         tbe = tb.TracebackException(
#             exc_type, exc_value, exc_tb,
#         )
#         log.error(logid+''.join(tbe.format()))


def write_unconstraint(
    save,
    sid,
    seq,
    unconstrained,
    data: PLFoldOutput,
    region,
    window,
    span,
    temperature,
    outdir,
    rawentry=None,
):
    """Fold Sequences with paired constraint

    Parameters
    ----------
    save: bool Should npy or npy and gz output be saved
    sid: str Id for Sequence
    seq: str Sequence to fold
    unconstrained: str Infix for ouput files of unconstrained folding
    data: PLFoldOutput Output of call to fold_unconstraint
    region: int -u Option of RNAplfold
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    temperature: float Temperature
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    rawentry: str Name for raw file to pass to write_unconstraint
    """
    logid = SCRIPTNAME + ".write_unconstraint: "
    log.debug(logid + " ".join([str(save), str(len(seq)), str(len(data.numpy_array)), str(rawentry)]))
    try:
        temperature = re.sub("[.,]", "", str(temperature))
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir, goi)

        try:
            gr = str(sid.split(":")[3].split("(")[0])
        except IndexError:
            gr = "na"

        if unconstrained != "STDOUT":
            if not os.path.exists(temp_outdir):
                try:
                    # Multiprocessing can lead to 'did not just yet exist but suddenly does' error and
                    # we do not want to catch that
                    os.makedirs(temp_outdir, exist_ok=True)
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
            if rawentry:
                filename = f"{goi}_{chrom}_{strand}_{rawentry}_{unconstrained}_{window}_{span}_{temperature}"
                gz_filepath = os.path.join(temp_outdir, f"{filename}.gz")
                npy_filepath = os.path.join(temp_outdir, f"{filename}.npy")
            else:
                filename = f"{goi}_{chrom}_{strand}_{gr}_{unconstrained}_{window}_{span}_{temperature}"
                gz_filepath = os.path.join(temp_outdir, f"{filename}.gz")
                npy_filepath = os.path.join(temp_outdir, f"{filename}.npy")
            if save > 0 and not os.path.exists(gz_filepath):
                with gzip.open(gz_filepath, "wb") as o:
                    out = data.get_text(nan="nan", truncated=True)
                    if out and len(out) > 1:
                        o.write(bytes(out, encoding="UTF-8"))
            if not os.path.exists(npy_filepath):
                printdiff(data.get_rnamediator_np_array(), npy_filepath)

        else:
            print(data.get_text(nan="nan", truncated=True))
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def write_constraint(
    save,
    sid,
    seq,
    paired,
    unpaired,
    data_u,
    data_p,
    constrain,
    region,
    diff_nu,
    diff_np,
    window,
    span,
    temperature,
    outdir,
):
    """Write output of constraint folding

    Parameters
    ----------
    save: bool Should npy or npy and gz output be saved
    sid: str Id for Sequence
    seq: str Sequence to fold
    paired: str Infix for output files of paired constraint
    unpaired: str Infix for ouput files of unpaired constraint
    data_u: PLFoldOutput Output of call to fold_unconstraint for unpaired constraint
    data_p: PLFoldOutput Output of call to fold_unconstraint for paired constraint
    constraint: str Id od constraint
    region: int -u Option of RNAplfold
    diff_nu: np.array Diff between unconstrained and constraint unpaired
    diff_np: np.array Diff between unconstrained and constraint paired
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    temperature: float Temperature
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    """
    logid = SCRIPTNAME + "write_constraint: "
    try:
        temperature = re.sub("[.,]", "", str(temperature))
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir, goi)
        # print outputs to file or STDERR
        if data_p:
            if paired != "STDOUT":
                if not os.path.exists(temp_outdir):
                    os.makedirs(temp_outdir)
                filename = f"StruCons_{goi}_{chrom}_{strand}_{constrain}_{paired}_{window}_{span}_{temperature}.gz"
                filepath = os.path.join(temp_outdir, filename)
                if save > 0 and not os.path.exists(filepath):
                    with gzip.open(filepath, "wb") as o:
                        out = data_p.get_text(nan="nan", truncated=True)
                        if out and len(out) > 1:
                            o.write(bytes(out, encoding="UTF-8"))
                        else:
                            log.error("No output produced " + sid)
            else:
                print(data_p.get_text(nan="nan", truncated=True))

        if unpaired != "STDOUT":
            if not os.path.exists(temp_outdir):
                os.makedirs(temp_outdir)
            filename = f"StruCons_{goi}_{chrom}_{strand}_{constrain}_{unpaired}_{window}_{span}_{temperature}.gz"
            filepath = os.path.join(temp_outdir, filename)
            if save > 0 and not os.path.exists(filepath):
                with gzip.open(filepath, "wb") as o:
                    out = data_u.get_text(nan="nan", truncated=True)
                    if out and len(out) > 1:
                        o.write(bytes(out, encoding="UTF-8"))
                    else:
                        log.error("No output produced " + sid)
        else:
            print(data_u.get_text(nan="nan", truncated=True))

        if not (diff_nu is None) and diff_nu.any():
            if unpaired != "STDOUT":
                if not os.path.exists(temp_outdir):
                    os.makedirs(temp_outdir)
                filename = f"StruCons_{goi}_{chrom}_{strand}_{constrain}_diffnu_{window}_{span}_{temperature}.npy"
                filepath = os.path.join(temp_outdir, filename)
                if not os.path.exists(filepath):
                    printdiff(diff_nu, filepath)
            else:
                _npprint(diff_nu)

        if not (diff_np is None) and diff_np.any():
            if unpaired != "STDOUT":
                if not os.path.exists(temp_outdir):
                    os.makedirs(temp_outdir)
                filename = f"StruCons_{goi}_{chrom}_{strand}_{constrain}_diffnp_{window}_{span}_{temperature}.npy"
                filepath = os.path.join(temp_outdir, filename)
                if not os.path.exists(filepath):
                    printdiff(diff_np, filepath)
            else:
                _npprint(diff_np)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


# def write_temp(save, sid, seq, temp, data, region, diff, window, span, outdir):
#
#     logid = SCRIPTNAME + '.write_temp: '
#     try:
#         logid = SCRIPTNAME + '.write_temp: '
#         #print outputs to file or STDERR
#         goi, chrom, strand = idfromfa(sid)
#         temp_outdir = os.path.join(outdir,goi)
#         if not os.path.exists(temp_outdir):
#             os.makedirs(temp_outdir)
#         if save > 0  and not os.path.exists(os.path.join(temp_outdir,'TempCons_'+goi+'_'+chrom+'_'+strand+'_'+temp+'_temp_'+window+'_'+str(span)+'.gz')):
#             with gzip.open(os.path.join(temp_outdir,'TempCons_'+goi+'_'+chrom+'_'+strand+'_'+temp+'_temp_'+window+'_'+str(span)+'.gz'), 'wb') as o:
#                 o.write(bytes(print_up(data['up'],len(seq),region),encoding='UTF-8'))
#             if not os.path.exists(os.path.join(temp_outdir,'TempCons_'+goi+'_'+chrom+'_'+strand+'_'+temp+'_difft_'+window+'_'+str(span)+'.gz')):
#                 with gzip.open(os.path.join(temp_outdir,'TempCons_'+goi+'_'+chrom+'_'+strand+'_'+temp+'_difft_'+window+'_'+str(span)+'.gz'), 'wb') as o:
#                     npprint(diff,o)
#         if not          os.path.exists(os.path.join(temp_outdir,'TempCons_'+goi+'_'+chrom+'_'+strand+'_'+temp+'_temp_'+window+'_'+str(span)+'.npy')):
#             printdiff(up_to_array(data['up']),os.path.join(temp_outdir,'TempCons_'+goi+'_'+chrom+'_'+strand+'_'+temp+'_temp_'+window+'_'+str(span)+'.npy'))
#
#     except Exception:
#         exc_type, exc_value, exc_tb = sys.exc_info()
#         tbe = tb.TracebackException(
#             exc_type, exc_value, exc_tb,
#             )
#         log.error(logid+''.join(tbe.format()))
#     return 1


def checkexisting(sid, paired, unpaired, cons, region, window, span, temperature, outdir):
    """Check if output file from constraint folding exists

    Parameters
    ----------
    sid: str Id for Sequence
    paired: str Infix for output files of paired constraint
    unpaired: str Infix for ouput files of unpaired constraint
    cons: str Id od constraint
    region: int -u Option of RNAplfold
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    temperature: float Temperature
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used

    Returns
    -------
    Bool: True/False if file exists/not
    """
    logid = SCRIPTNAME + ".checkexisting: "
    try:
        goi, chrom, strand = idfromfa(sid)
        log.debug(f"{logid} ID:{sid}, GOI:{goi}, OUT:{outdir}")
        temp_outdir = os.path.join(outdir, goi)
        temperature = re.sub("[.,]", "", str(temperature))
        p = os.path.join(
            temp_outdir,
            "StruCons_"
            + goi
            + "_"
            + chrom
            + "_"
            + strand
            + "_"
            + cons
            + "_"
            + str(paired)
            + "_"
            + str(window)
            + "_"
            + str(span)
            + "_"
            + str(temperature)
            + ".gz",
        )
        u = os.path.join(
            temp_outdir,
            "StruCons_"
            + goi
            + "_"
            + chrom
            + "_"
            + strand
            + "_"
            + cons
            + "_"
            + str(unpaired)
            + "_"
            + str(window)
            + "_"
            + str(span)
            + "_"
            + str(temperature)
            + ".gz",
        )
        log.debug(f"{logid}: checking for {p} and {u}")

        if paired:
            if os.path.exists(p) and os.path.exists(u):
                return True
        elif os.path.exists(u):
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


def check_raw_existing(sid, unconstrained, cons, region, window, span, temperature, outdir, rawentry=None):
    """Check if output file from unconstrained folding exists

    Parameters
    ----------
    sid: str Id for Sequence
    unconstrained: str Infix for ouput files of unconstrained folding
    cons: str Id od constraint
    region: int -u Option of RNAplfold
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    temperature: float Temperature
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    rawentry: str Name for raw file to pass to write_unconstraint

    Returns
    -------
    Bool: True/False if file exists/not
    """
    logid = SCRIPTNAME + ".check_raw_existing: "
    try:
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir, goi)
        temperature = re.sub("[.,]", "", str(temperature))
        try:
            gr = str(sid.split(":")[3].split("(")[0])
        except IndexError:
            gr = "na"
        if rawentry:
            filename = f"{goi}_{chrom}_{strand}_{rawentry}_{unconstrained}_{window}_{span}_{temperature}"
            gz_filepath = os.path.join(temp_outdir, f"{filename}.gz")
            npy_filepath = os.path.join(temp_outdir, f"{filename}.npy")
        else:
            filename = f"{goi}_{chrom}_{strand}_{gr}_{unconstrained}_{window}_{span}_{temperature}"
            gz_filepath = os.path.join(temp_outdir, f"{filename}.gz")
            npy_filepath = os.path.join(temp_outdir, f"{filename}.npy")
        log.debug(f"{logid} Checking {os.path.join(temp_outdir, npy_filepath)}")
        if os.path.exists(os.path.join(temp_outdir, npy_filepath)):
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


def main(args=None):
    """Main process, prepares run_settings dict, creates logging process queue and worker processes for folding, calls fold

    Parameters
    ----------
    args: str, optional
        Arguments for main run

    Returns
    -------
    Call to pl_fold
    """

    logid = SCRIPTNAME + ".main: "
    try:
        if not args:
            args = parseargs_plcons()

        log.debug(f"{args}")

        if args.version:
            sys.exit("Running RNAmediator version " + __version__)

        queue, listener, worker_configurer = rnamediator_logging_setup(args.logdir, args.loglevel, SCRIPTNAME)

        log.info(logid + "Running " + SCRIPTNAME + " on " + str(args.procs) + " cores.")

        log.info(logid + "CLI: " + sys.argv[0] + " " + "{}".format(" ".join([shlex.quote(s) for s in sys.argv[1:]])))

        run_settings, outdir = preprocess(
            args.sequence, args.constrain, args.conslength, args.constype, args.outdir, args.genes
        )

        pl_fold(
            args.window,
            args.span,
            args.region,
            args.temperature,
            args.multi,
            args.unconstrained,
            args.unpaired,
            args.paired,
            args.save,
            args.procs,
            outdir,
            run_settings,
            queue=queue,
            configurer=worker_configurer,
            level=args.loglevel,
        )
        queue.put(None)
        queue.join()
        listener.join()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        if log:
            log.error(logid + "".join(tbe.format()))
        else:
            print(f'ERROR: {logid} {"".join(tbe.format())}')
        queue.terminate()
        listener.join()
        sys.exit(1)


####################
####    MAIN    ####
####################

if __name__ == "__main__":

    outer_logid = SCRIPTNAME + ".main: "
    try:
        main()

    except Exception:
        outer_exc_type, outer_exc_value, outer_exc_tb = sys.exc_info()
        outer_tbe = tb.TracebackException(
            outer_exc_type,
            outer_exc_value,
            outer_exc_tb,
        )
        if log:
            log.error(outer_logid + "".join(outer_tbe.format()))
        else:
            print(f'ERROR: {logid} {"".join(outer_tbe.format())}')

    # ConstraintPLFold.py ends here
