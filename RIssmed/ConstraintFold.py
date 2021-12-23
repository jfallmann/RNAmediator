#!/usr/bin/env python3
### ConstraintFold.py ---
##
## Filename: ConstraintFold.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Sep  6 09:02:18 2018 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Mon Jan 25 10:23:45 2021 (+0100)
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
# other modules
import os
import sys
import importlib
import multiprocessing
import shlex
import errno

# ViennaRNA
import RNA

from RIssmed import _version

__version__ = _version.get_versions()["version"]

# load own modules
from RIssmed.Tweaks.FileProcessor import *
from RIssmed.Tweaks.RNAtweaks import (
    _get_bppm as get_bppm,
    _get_ddg as get_ddg,
    _calc_gibbs as calc_gibbs,
    _calc_ddg as calc_ddg,
    _calc_bpp as calc_bpp,
    _calc_nrg as calc_nrg,
    _constrain_paired as constrain_paired,
    _constrain_unpaired as constrain_unpaired,
    get_location,
    api_rnafold,
    FoldOutput,
)
from RIssmed.Tweaks.NPtweaks import *
from RIssmed.Tweaks.RIssmed import (
    preprocess,
    SequenceSettings,
    rissmed_logging_setup,
    expand_window,
    localize_window,
)

log = logging.getLogger(__name__)  # use module name
SCRIPTNAME = os.path.basename(__file__).replace(".py", "")


def fold(
    run_settings,
    outdir,
    window,
    span,
    temp,
    constrain,
    conslength,
    procs,
    save="STDOUT",
    pattern=None,
    cutoff=None,
    queue=None,
    configurer=None,
    level=None,
):
    """fold prepares and submits sequences to be folded constrained and unconstrained

    Parameters
    ----------
    run_settings: Dict[str, SequenceSettings] RIssmed run settings dictionary using fasta ids as keys and Sequence Settings as values
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    temp: int
        Temperature to fold at
    constrain : str The file location of constrain file
    conslength : int Length of the constraint, only used if constrain is sliding
    procs: int Number of processes to run in parallel
    save: str
        The name of the output file to generate or STDOUT
    pattern: str String pattern for gene of interest
    cutoff: float Cutoff for raw accessibility, regions below this propability of being unpaired will not be folded
    queue: multiprocessing_queue Logging process queue
    configurer: multiprocessing_config for Logging processes
    level: logging.level Level for log process

    Returns
    -------
    Call to constrain_seq
    """

    logid = SCRIPTNAME + ".fold: "
    try:
        if queue and level:
            configurer(queue, level)

        # Create process pool with processes
        num_processes = procs or 1
        # with get_context("spawn").Pool(processes=num_processes-1, maxtasksperchild=1) as pool:
        pool = multiprocessing.Pool(processes=num_processes, maxtasksperchild=1)

        # Start the work
        for fasta in run_settings:
            fasta_settings = run_settings[fasta]
            log.debug(logid + str(dir(fasta_settings)))

            goi = fasta_settings.gene
            gs = fasta_settings.genomic_coords.start
            ge = fasta_settings.genomic_coords.end
            gstrand = fasta_settings.genomic_coords.strand
            seq_record = fasta_settings.sequence_record
            conslist = fasta_settings.constrainlist

            log.debug(
                logid
                + "Constraints: "
                + str(conslist)
                + " Seq_record: "
                + str(seq_record.seq)
            )
            fa = str(seq_record.seq)

            if len(fa) < window:
                log.warning(
                    str(
                        "Sequence of "
                        + goi
                        + " too short, seqlenght "
                        + str(len(fa))
                        + " with window size "
                        + str(window)
                    )
                )
                continue

            if pattern and pattern not in goi:
                continue

            log.info(logid + "Working on " + goi + "\t" + seq_record.id)

            # constraints
            for cons_tuple in conslist:
                log.debug(logid + "ENTRY: " + str(cons_tuple))
                cons_tuple = [str(cons) for cons in cons_tuple]
                cons = ":".join(cons_tuple)

                if not window:
                    window = len(fa)
                if not span:
                    span = len(fa)

                if (
                    cons == "NOCONS"
                ):  # in case we just want to fold the sequence without constraints at all
                    gibbs_uc = [
                        pool.apply_async(
                            fold_unconstraint,
                            args=(fa),
                            kwds={
                                "temp": temp,
                                "queue": queue,
                                "configurer": configurer,
                                "level": level,
                            },
                        )
                    ]
                    return gibbs_uc

                else:
                    # we now have a list of constraints and for the raw seq comparison we only need to fold windows around these constraints
                    fstart, fend = [None, None]
                    start, end = [None, None]

                    if len(cons_tuple) > 1 or (
                        any(x in cons for x in ["paired", "Paired"]) or ":" in cons
                    ):  # Not strand dependend, still genomic coords; Paired constraints should be used

                        if gstrand == "+" or gstrand == ".":
                            [fstart, fend], [start, end] = [
                                [x - gs for x in get_location(cn)[:2]]
                                for cn in cons.split(":", 1)
                            ]
                        else:
                            [fstart, fend], [start, end] = [
                                [ge - x for x in get_location(cn)[:2][::-1]]
                                for cn in cons.split(":", 1)
                            ]
                        cons = (
                            str(fstart)
                            + "-"
                            + str(fend)
                            + ":"
                            + str(start)
                            + "-"
                            + str(end)
                        )
                        consstr = str(f"{fstart}-{fend}:{start}-{end}")
                        if start < 0 or fstart < 0 or end > len(fa) or fend > len(fa):
                            log.warning(
                                logid
                                + "Constraint out of sequence bounds! skipping! "
                                + ",".join(
                                    map(
                                        str,
                                        [
                                            goi,
                                            len(fa),
                                            str(fstart) + "-" + str(fend),
                                            str(start) + "-" + str(end),
                                        ],
                                    )
                                )
                            )
                            continue

                        # if checkexisting(
                        #    str(seq_record.id), paired, unpaired, cons, region, window, span, outdir
                        # ):
                        #    log.warning(
                        #        logid + str(cons) + ' Exists for ' + str(seq_record.id) + '! #Skipping!'
                        #    )
                        #    continue
                    else:
                        if gstrand == "+" or gstrand == ".":
                            fstart, fend = [x - gs for x in get_location(cons)[:2]]
                        else:
                            fstart, fend = [
                                ge - x for x in get_location(cons)[:2][::-1]
                            ]

                        consstr = str(fstart) + "-" + str(fend)

                        log.debug(logid + str.join(" ", [goi, consstr, gstrand]))

                    genecoords = list([gs, ge, gstrand])
                    const = list([fstart, fend, start, end])

                    pool.apply_async(
                        constrain_seq,
                        args=(
                            seq_record,
                            const,
                            conslength,
                            consstr,
                            window,
                            span,
                            temp,
                            save,
                            outdir,
                            genecoords,
                        ),
                        kwds={"queue": queue, "configurer": configurer, "level": level},
                    )

            pool.close()
            pool.join()

        log.info(logid + "DONE: output in: " + str(outdir))

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


##### Functions #####


def constrain_seq(
    seq_record,
    const,
    conslength,
    cons,
    window,
    span,
    temp,
    save="STDOUT",
    outdir=None,
    genecoords=None,
    queue=None,
    configurer=None,
    level=None,
):
    """takes input of fold and folds constrained and unconstrained sequences

    Parameters
    ----------
    seq_record : SeqIO.SeqRecord Sequence record of the sequence to which the constraint should be added
    const : list The constraint position(s) as list
    conslength : int Length of the constraint, only used if constrain is sliding
    cons : str The constraint annotation
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    temp: int
        Temperature to fold at
    save: str
        Name of output file or STDOUT
    outdir : str Location for the output directory
    genecoords: list Genomic coordinates of gene of interest
    procs: int Number of processes to run in parallel
    queue: multiprocessing_queue Logging process queue
    configurer: multiprocessing_config for Logging processes
    level: logging.level Level for log process

    Returns
    -------
    Call to writeout which prints results to file
    """

    logid = SCRIPTNAME + ".constrain_seq: "
    try:
        if queue and level:
            configurer(queue, level)

        seq = str(seq_record.seq).upper().replace("T", "U")
        goi, chrom, strand = idfromfa(seq_record.id)
        log.debug(logid + " Sequence to fold: " + str(seq))

        if not window:
            window = len(seq)
        if not span:
            span = window

        if len(seq) < int(window):
            log.error(
                logid + "Sequence to small, skipping " + str(id) + "\t" + str(len(seq))
            )
            return

        if genecoords and goi and goi in genecoords:
            gs, ge, gstrand = get_location(genecoords[goi][0])
            if gstrand != strand:
                log.warning(
                    logid
                    + "Strand values differ between Gene annotation and FASTA file! Please check your input for "
                    + str(goi)
                )
        else:
            gs = ge = 0
            gstrand = "."
            log.warning(
                logid
                + "No coords found for gene "
                + goi
                + "! Assuming coordinates are already local!"
            )

        dist = None
        fstart, fend, start, end = const

        if end:
            dist = abs(fstart - end)
            if dist > window:
                return log.warning(
                    logid
                    + "Window "
                    + str(window)
                    + " too small for constraint distance "
                    + str(dist)
                )

        # we no longer fold the whole sequence but only the constraint region +- window size
        if window < len(seq):
            if dist:
                tostart, toend = expand_window(fstart, end, window, len(seq), dist)
                log.debug(
                    logid
                    + "Considering distance "
                    + str(dist)
                    + " between constraints "
                    + str(const)
                    + " for window extraction from "
                    + str(tostart)
                    + " to "
                    + str(toend)
                )
            else:
                tostart, toend = expand_window(fstart, fend, window, len(seq))

            if fend is not None and toend < fend:
                log.warning(
                    logid
                    + "Constraint "
                    + str(cons)
                    + " out of sequence range "
                    + str(toend)
                    + "!Skipping!"
                )  # One of the constraints is outside the sequence window
                return
        else:
            tostart, toend = [1, len(seq)]

        seqtofold = str(seq[tostart - 1 : toend])

        log.debug(
            logid
            + ",".join(
                [
                    str(start),
                    str(end),
                    str(tostart),
                    str(toend),
                    "SEQUENCE: " + seqtofold,
                ]
            )
        )

        if len(seqtofold) < window:
            log.warning(
                logid
                + "Sequence of "
                + goi
                + " to short, seqlenght "
                + str(len(seqtofold))
                + " with window size "
                + str(window)
                + "!Skipping! "
            )
            raise Exception(
                "Sequence of "
                + goi
                + " to short, seqlenght "
                + str(len(seqtofold))
                + " with window size "
                + str(window)
            )

        # get local start,ends 0 based closed
        cfstart, cfend = localize_window(fstart, fend, tostart, toend, len(seq))
        check = []

        if start and end:
            cstart = start - tostart
            cend = end - tostart
            check = (cfstart, cfend, cstart, cend)
        else:
            check = (cfstart, cfend)

        # Coordinates for fold annotation
        coords = [gs, ge, tostart, toend]

        sp = ep = None
        os, oe = [fstart + tostart + gs - 1, fend + tostart + gs]
        gtostart, gtoend = [tostart + gs, toend + gs]
        printcons = str.join(
            "|",
            [
                str.join("-", [str(tostart), str(toend)]),
                str.join("-", [str(gtostart), str(gtoend)]),
                str.join("-", [str(fstart + 1), str(fend + 1)]),
            ],
        )

        if len(check) > 2:
            osp, oep = [start + tostart + gs - 1, end + tostart + gs]
            printcons = printcons + ":" + str.join("-", [str(start), str(end)])

        printcons = printcons + "|" + str.join("-", [str(os), str(oe)])
        if len(check) > 2:
            printcons = printcons + ":" + str.join("-", [str(osp), str(oep)])

        Output = FoldOutput()
        Output = api_rnafold(
            seqtofold, span, 37, None, FoldOut=Output, coordinates=coords
        )

        for cons in ["paired", "unpaired"]:
            checkc = (cons,) + check
            Output = api_rnafold(
                seqtofold,
                span,
                37,
                constraint=[checkc],
                FoldOut=Output,
                coordinates=coords,
            )

        if len(check) > 3:
            for cons in [
                "secondconstraint_paired",
                "secondconstraint_unpaired",
                "bothconstraint_paired",
                "bothconstraint_unpaired",
            ]:
                checkc = (cons,) + check
                Output = api_rnafold(
                    seqtofold,
                    span,
                    37,
                    constraint=[checkc],
                    FoldOut=Output,
                    coordinates=coords,
                )

        Output.annotate(
            str(goi)
            + ": "
            + str.join(",", [str(chrom), str(gs), str(ge), str(gstrand)])
        )

        write_out(Output, window, span, printcons, save, outdir)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def constrain_temp(
    fa, temp, window, span, an, save, outdir, queue=None, configurer=None, level=None
):
    """Takes Sequence and fold under temperature constraint
    Currently not implemented

    Parameters
    ----------
    fa: SeqIO.SeqRecord Sequence record of the sequence to which the constraint should be added
    temp : int Temperature[K]
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    an: str Name for outfile
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    queue: multiprocessing_queue Logging process queue
    configurer: multiprocessing_config for Logging processes
    level: logging.level Level for log process

    Returns
    -------
    NotImplementedError()
    """

    # refresh model details
    raise NotImplementedError("Needs to be reimplemented")

    logid = SCRIPTNAME + ".constrain_temp: "
    log.info(logid + "Constraining Temp to " + temp)
    try:
        if queue and level:
            configurer(queue, level)

        md = RNA.md()
        md.max_bp_span = span
        md.window_size = window
        # set temperature
        md.temperature = temp
        # create new fold_compound objects
        fc_t = RNA.fold_compound(str(seq_record.seq), md, RNA.OPTION_WINDOW)
        # data
        data = {"seq": str(seq_record.seq), "stru": "", "nrg": ""}
        # set model details
        md = RNA.md()
        # create new fold_compound object
        fc = RNA.fold_compound(data["seq"], md, RNA.OPTION_PF)
        gibbs = fc.pf()
        bppm = fc.bpp()

        data["stru"] = gibbs[0]
        data["nrg"] = gibbs[1]

        # call prop window calculation
        log.debug(logid + ";".join(map(str, [gibbs, constrain, conslength])))

        for item in bppm:
            for i in range(int(constrain), int(constrain) + conslength):
                log.debug(bppm.index(item), i, item[i])

        at = up_to_array(data_t["up"], None, len(seq_record.seq))

        diff_nt = an - at

        #####Set temp, and rewrite save for temp
        if save:
            log.info(logid + "SAVINGTEMP")
            write_temp(fa, str(temp), data_t, diff_nt, str(window), outdir)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def foldaround(
    seq, fc, pos, clength, gibbs, nrg, queue=None, configurer=None, level=None
):
    """here we take the sequence and the RNA.fold_compound and constrain regions of length clength around it to see what happens at the original binding site

    Parameters
    ----------
    seq: SeqIO.SeqRecord Sequence record of the sequence to which the constraint should be added
    fc : RNA.fold_compound
    pos : int constraint start
    clength : int The constraint length
    gibbs: float ddG
    nrg: float Folding energy
    queue: multiprocessing_queue Logging process queue
    configurer: multiprocessing_config for Logging processes
    level: logging.level Level for log process

    Returns
    -------
    FoldOutput object
    """

    logid = SCRIPTNAME + ".foldaround: "
    try:
        if queue and level:
            configurer(queue, level)

        cstart = pos
        cend = pos + clength - 1
        fc = constrain_unpaired(fc, cstart, cend)
        gibbs_u = calc_gibbs(fc)
        bppm = get_bppm(fc.bpp(), cstart, cend)
        bpp = calc_bpp(bppm)
        ddg = gibbs_u - gibbs
        nrg_u = calc_nrg(bpp)
        nrg_diff = nrg_u - nrg

        return [gibbs_u, ddg, nrg_diff]

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def fold_unconstraint(seq, queue=None, configurer=None, level=None):
    """here we take the sequence and the RNA.fold_compound and fold without constraint

    Parameters
    ----------
    seq: SeqIO.SeqRecord Sequence record of the sequence to which the constraint should be added
    queue: multiprocessing_queue Logging process queue
    configurer: multiprocessing_config for Logging processes
    level: logging.level Level for log process

    Returns
    -------
    FoldOutput object
    """

    logid = SCRIPTNAME + ".fold_unconstraint: "
    try:
        if queue and level:
            configurer(queue, level)

        fold_output = api_rnafold(seq)

        return fold_output
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def write_out(Output, window, span, const, fname="STDOUT", outdir=None):
    """Print results to file

    Parameters
    ----------
    Output: FoldOutput object
        Contains all results generated by constrain
    window: int
        Window to fold
    span: int
        Max basepair span
    cons: str
        Constraint applied on fold
    fname: str, optional
        Name of file to write to or STDOUT
    outdir: str, optional
        Name of output directory or None

    Returns
    -------
    Print results to file or STDOUT
    """

    try:
        logid = SCRIPTNAME + ".write_out: "
        goi, chrom, start, end, strand = Output.parse_anno()
        temp_outdir = os.path.join(outdir, goi)

        if fname != "STDOUT":
            if not os.path.exists(temp_outdir):
                os.makedirs(temp_outdir)
            if not os.path.exists(
                os.path.join(
                    temp_outdir,
                    "_".join([goi, chrom, strand, fname, const, str(window), str(span)])
                    + ".gz",
                )
            ):
                with gzip.open(
                    os.path.join(
                        temp_outdir,
                        "_".join(
                            [goi, chrom, strand, fname, const, str(window), str(span)]
                        )
                        + ".gz",
                    ),
                    "wb",
                ) as o:
                    o.write(bytes(Output.get_text() + "\n", encoding="UTF-8"))
            else:
                log.info(
                    os.path.join(
                        temp_outdir,
                        "_".join(
                            [goi, chrom, strand, fname, const, str(window), str(span)]
                        )
                        + ".gz",
                    )
                    + " exists, will append!"
                )
                with gzip.open(
                    os.path.join(
                        temp_outdir,
                        "_".join(
                            [goi, chrom, strand, fname, const, str(window), str(span)]
                        )
                        + ".gz",
                    ),
                    "ab",
                ) as o:
                    o.write(bytes(Output.get_text(h=False), encoding="UTF-8"))
        else:
            print(Output.get_text() + "\n")

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def checkexisting(
    sid,
    paired,
    unpaired,
    cons,
    window,
    span,
    outdir,
    queue=None,
    configurer=None,
    level=None,
):
    """Validates if output file already exists

    Parameters
    ----------
    sid: SeqIO.SeqID Sequence ID of the sequence to which the constraint should be added
    unpaired: str Suffix for ouput files of unpaired constraint
    paired: str Suffix for output files of paired constraint
    cons : str The annoation for the constraint
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    queue: multiprocessing_queue Logging process queue
    configurer: multiprocessing_config for Logging processes
    level: logging.level Level for log process

    Returns
    -------
    True/False if file exists
    """

    logid = SCRIPTNAME + ".checkexisting: "
    try:
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir, goi)

        if os.path.exists(
            os.path.join(
                temp_outdir,
                "_".join([goi, chrom, strand, fname, const, window, span]) + ".gz",
            )
        ):
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
    return 1


def main(args=None):
    """Main process, prepares run_settings dict, creates logging process queue and worker processes for folding, calls fold

    Parameters
    ----------

    Returns
    -------
    Call to fold
    """

    logid = SCRIPTNAME + ".main: "
    try:
        if not args:
            args = parseargs_foldcons()

        if args.version:
            sys.exit("Running RIssmed version " + __version__)

        queue, listener, worker_configurer = rissmed_logging_setup(
            args.logdir, args.loglevel, SCRIPTNAME
        )

        log.info(logid + "Running " + SCRIPTNAME + " on " + str(args.procs) + " cores.")
        log.info(
            logid
            + "CLI: "
            + sys.argv[0]
            + " "
            + "{}".format(" ".join([shlex.quote(s) for s in sys.argv[1:]]))
        )

        run_settings, outdir = preprocess(
            args.sequence, args.constrain, args.conslength, args.outdir, args.genes
        )

        fold(
            run_settings,
            outdir,
            args.window,
            args.span,
            args.temperature,
            args.constrain,
            args.conslength,
            args.procs,
            args.save,
            args.pattern,
            args.cutoff,
            queue=queue,
            configurer=worker_configurer,
            level=args.loglevel,
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
    SCRIPTNAME = os.path.basename(__file__).replace(".py", "")
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

# ConstraintFold.py ends here
