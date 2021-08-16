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
import importlib
import multiprocessing
import shlex

# Biopython stuff
from Bio import SeqIO

# numpy
# RNA
# Logging
import datetime
from Tweaks.logger import (
    makelogdir,
    makelogfile,
    listener_process,
    listener_configurer,
    worker_configurer,
)

# load own modules
from Tweaks.FileProcessor import *
from Tweaks.RNAtweaks import *
from Tweaks.NPtweaks import *
from Tweaks.RIssmed import (
    preprocess,
    SequenceSettings,
    rissmed_logging_setup,
    expand_window,
    localize_window,
)

log = logging.getLogger(__name__)  # use module name
SCRIPTNAME = os.path.basename(__file__).replace('.py', '')


def fold(
    run_settings,
    outdir,
    window,
    span,
    unconstraint,
    unpaired,
    paired,
    constrain,
    conslength,
    save,
    procs,
    pattern=None,
    cutoff=None,
    queue=None,
    configurer=None,
    level=None,
):
    """fold prepares and submitts sequenced to be folded constrained and unconstrained

    Parameters
    ----------
    run_settings: Dict[str, SequenceSettings] RIssmed run settings dictionary using fasta ids as keys and Sequence Settings as values
    outdir : str Location of the Outpu directory. If it is an empty string os.cwd() is used
    window: int Size of window to fold
    span: int Maximum base-pair span to be evaluated
    unconstraint: str Suffix for unconstraint output files
    unpaired: str Suffix for ouput files of unpaired constraint
    paired: str Suffix for output files of paired constraint
    constrain : str The file location of constrain file
    conslength : int Length of the constraint, only used if constrain is sliding
    save: bool If output files should be saved
    procs: int Number of processes to run in parallel
    pattern: str String pattern for gene of interest
    cutoff: float Cutoff for raw accessibility, regions below this propability of being unpaired will not be folded
    queue: multiprocessing_queue Logging process queue
    configurer: multiprocessing_config for Logging processes
    level: logging.level Level for log process

    Returns
    -------
    Call to constrain_seq
    """

    logid = SCRIPTNAME + '.fold: '
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
            goi = fasta_settings.gene
            gs = fasta_settings.genomic_coords.start
            ge = fasta_settings.genomic_coords.end
            gstrand = fasta_settings.genomic_coords.strand
            fa = fasta_settings.sequence_record
            conslist = fasta_settings.constrainlist
            log.debug(logid + str(conslist))

            if len(seq_record.seq) < window:
                log.warning(
                    str(
                        'Sequence of '
                        + goi
                        + ' too short, seqlenght '
                        + str(len(seq_record.seq))
                        + ' with window size '
                        + str(window)
                    )
                )
                continue

            if pattern and pattern not in goi:
                continue

            log.info(logid + 'Working on ' + goi + "\t" + seq_record.id)

            # constraints
            for cons_tuple in conslist:
                log.debug(logid + 'ENTRY: ' + str(cons_tuple))

                if not window:
                    window = len(seq_record.seq)
                if not span:
                    span = len(seq_record.seq)

                if entry == 'NOCONS':  # in case we just want to fold the sequence without constraints at all
                    gibbs_uc = [
                        pool.apply_async(
                            fold_unconstraint,
                            args=(seq_record.seq),
                            kwds={'queue': queue, 'configurer': configurer, 'level': level},
                        )
                    ]
                    return gibbs_uc

                else:
                    # we now have a list of constraints and for the raw seq comparison we only need to fold windows around these constraints
                    fstart, fend = [None, None]
                    start, end = [None, None]

                    if len(cons_tuple) > 1 or (
                        any(x in constrain for x in ['paired', 'Paired']) or ':' in cons_tuple
                    ):  # Not strand dependend, still genomic coords; Paired constraints should be used
                        cons_tuple = [str(cons) for cons in cons_tuple]
                        cons = ":".join(cons_tuple)

                        if gstrand == '+' or gstrand == '.':
                            [fstart, fend], [start, end] = [
                                [x - gs for x in get_location(cn)[:2]] for cn in cons.split(':', 1)
                            ]
                        else:
                            [fstart, fend], [start, end] = [
                                [ge - x for x in get_location(cn)[:2][::-1]] for cn in cons.split(':', 1)
                            ]
                        cons = str(fstart) + '-' + str(fend) + ':' + str(start) + '-' + str(end)
                        if start < 0 or fstart < 0 or end > len(seq_record.seq) or fend > len(seq_record.seq):
                            log.warning(
                                logid
                                + 'Constraint out of sequence bounds! skipping! '
                                + ','.join(
                                    map(
                                        str,
                                        [
                                            goi,
                                            len(seq_record.seq),
                                            str(start) + '-' + str(end),
                                            str(fstart) + '-' + str(fend),
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
                        if gstrand == '+' or gstrand == '.':
                            start, end = [x - gs for x in get_location(entry)[:2]]
                        else:
                            start, end = [ge - x for x in get_location(entry)[:2][::-1]]

                        cons = str(start) + '-' + str(end)

                        log.debug(logid + str.join(' ', [goi, cons, gstrand]))

                    genecoords = list(gs, ge, gstrand)
                    const = list(start, end, fstart, fend)
                    pool.apply_async(
                        constrain_seq,
                        args=(
                            seq_record,
                            const,
                            conslength,
                            cons,
                            window,
                            span,
                            unconstraint,
                            paired,
                            unpaired,
                            save,
                            outdir,
                            genecoords,
                        ),
                        kwds={'queue': queue, 'configurer': configurer, 'level': level},
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
        log.error(logid + ''.join(tbe.format()))


##### Functions #####


def constrain_seq(
    seq_record,
    const,
    conslength,
    cons,
    window,
    span,
    unconstraint,
    paired,
    unpaired,
    save,
    outdir,
    genecoords,
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
    unconstraint: str Suffix for unconstraint output files
    unpaired: str Suffix for ouput files of unpaired constraint
    paired: str Suffix for output files of paired constraint
    save: bool If output files should be saved
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

    logid = SCRIPTNAME + '.constrain_seq: '
    try:
        if queue and level:
            configurer(queue, level)

        outlist = list()
        goi, chrom, strand = idfromfa(seq_record.id)
        seq = seq_record.seq.upper().replace("T", "U")

        if not window:
            window = len(seq)
        if not span:
            span = window

        dist = None
        if fend:
            dist = abs(start - fend)
            if dist > window:
                return log.warning(
                    logid + 'Window ' + str(window) + ' too small for constraint distance ' + str(dist)
                )

        start, end, fstart, fend = const

        # we no longer fold the whole sequence but only the constraint region +- window size
        if window is not None:
            if dist:
                tostart = start - (window - abs(dist))
                toend = fend + (window - abs(dist)) + 1
                log.debug(
                    logid
                    + 'Considering distance '
                    + str(dist)
                    + ' between constraints '
                    + str(const)
                    + ' for window extraction from '
                    + str(tostart)
                    + ' to '
                    + str(toend)
                )
            else:
                tostart = start - window
                toend = end + window + 1

            if tostart < 0:
                tostart = 0
            if toend > len(seq):
                toend = len(seq)
            if fend is not None and toend < fend:
                log.warning(
                    logid + 'Constraint ' + str(cons) + ' out of sequence range ' + str(toend) + '!Skipping!'
                )  # One of the constraints is outside the sequence window
                return
            seqtofold = str(seq[tostart - 1 : toend])

        log.debug(
            logid + ','.join([str(start), str(end), str(tostart), str(toend), 'SEQUENCE: ' + seqtofold])
        )

        if len(seqtofold) < window:
            log.warning(
                logid
                + 'Sequence of '
                + goi
                + ' to short, seqlenght '
                + str(len(seqtofold))
                + ' with window size '
                + str(window)
                + '!Skipping! '
            )
            raise Exception(
                'Sequence of '
                + goi
                + ' to short, seqlenght '
                + str(len(seqtofold))
                + ' with window size '
                + str(window)
            )

        cstart = start - tostart
        cend = end - tostart

        if cstart < 0 or cend > len(seq):
            log.warning(
                logid
                + 'start of constraint '
                + str(cstart)
                + ' end of constraint '
                + str(cend)
                + ' while length of sequence '
                + str(len(seq))
                + '! Skipping!'
            )
            return

        checklist = []

        checklist.append((cstart, cend))
        if fstart and fend:
            cfstart = fstart - tostart
            cfend = fend - tostart
            checklist.append((cfstart, cfend))
            checklist.append((cstart, cend, cfstart, cfend))

        # data
        data = {'seq': seqtofold, 'stru': []}

        # set model details
        md = RNA.md()
        md.max_bp_span = span

        log.debug(logid + 'Constraints for ' + goi + ' are ' + str(checklist))

        for check in checklist:
            s, e, os, oe = [None, None, None, None]
            if genecoords:
                if goi in genecoords:
                    gs, ge, gstrand = get_location(genecoords[goi][0])
                    if gstrand != strand:
                        log.warning(
                            logid
                            + 'Strand values differ between Gene annotation and FASTA file! Please check your input for '
                            + str(goi)
                        )
                else:
                    gs, ge, gstrand = 0, 0, '.'
                    log.warning(
                        logid
                        + 'No coords found for gene '
                        + goi
                        + '! Assuming coordinates are already local!'
                    )
            else:
                gs = ge = 0
                gstrand = '.'
                log.warning(
                    logid + 'No coords found for gene ' + goi + '! Assuming coordinates are already local!'
                )

            log.debug(
                logid + 'GENECOORDS: ' + str(goi) + ': ' + str.join(',', [str(gs), str(ge), str(gstrand)])
            )

            if len(check) < 3:
                s, e = check
                sp, ep = [s + tostart + gs - 1, e + tostart + gs]
                gtostart, gtoend = [tostart + gs, toend + gs]
                printcons = str.join(
                    '|',
                    [
                        str.join('-', [str(tostart), str(toend)]),
                        str.join('-', [str(gtostart), str(gtoend)]),
                        str.join('-', [str(s), str(e)]),
                        str.join('-', [str(sp), str(ep)]),
                    ],
                )
            else:
                s, e, os, oe = check
                sp, ep, osp, oep = [
                    s + tostart + gs - 1,
                    e + tostart + gs,
                    os + tostart + gs - 1,
                    oe + tostart + gs,
                ]
                gtostart, gtoend = [tostart + gs, toend + gs]
                log.debug(
                    logid
                    + 'PAIRED:'
                    + ';'.join(
                        map(str, [s, e, os, oe, gs, ge, sp, ep, osp, oep, tostart, toend, gtostart, gtoend])
                    )
                )
                printcons = str.join(
                    '|',
                    [
                        str.join('-', [str(tostart), str(toend)]),
                        str.join('-', [str(gtostart), str(gtoend)]),
                        str.join(':', [str.join('-', [str(s), str(e)]), str.join('-', [str(os), str(oe)])]),
                        str.join(
                            ':', [str.join('-', [str(sp), str(ep)]), str.join('-', [str(osp), str(oep)])]
                        ),
                    ],
                )

            # create new fold_compound object
            fc = RNA.fold_compound(data['seq'], md)
            fc_u = RNA.fold_compound(data['seq'], md)
            fc_p = RNA.fold_compound(data['seq'], md)

            # call pf and prop calculation
            gibbs = fc.pf()[1]
            bppm = get_bppm(fc.bpp(), cstart, cend)

            if bppm is None:
                log.error(logid + 'Empty bpp matrix returned, stopping here!')
                #                sys.exit(logid+'Empty bpp matrix returned, stopping here!')
                return

            # enforce paired
            fc_p = constrain_paired(fc_p, s, e)
            # enforce unpaired
            fc_u = constrain_unpaired(fc_u, s, e)
            # calculate probs and nrg
            gibbs_u = calc_gibbs(fc_u)
            bppm_u = get_bppm(fc_u.bpp(), cstart, cend)
            dg_u = gibbs_u - gibbs

            gibbs_p = calc_gibbs(fc_p)
            bppm_p = get_bppm(fc_p.bpp(), cstart, cend)
            dg_p = gibbs_p - gibbs

            ###  access_energy([a, b]) = -RT log(prob([a, b]))
            ###  prob([a, b]) = sum_{s in S[a, b]} exp(-E(s)/RT) / sum_{s in S0} exp(-E(s)/RT)
            bpp = calc_bpp(bppm)
            bpp_u = calc_bpp(bppm_u)
            bpp_p = calc_bpp(bppm_p)
            nrg = calc_nrg(bpp)
            nrg_u = calc_nrg(bpp_u)
            nrg_p = calc_nrg(bpp_p)

            fn = 'constraint'

            if os and oe:
                fn = 'pairedconstraint'

            outlist.append(
                [fa, fn, gibbs, '0', nrg, printcons, str(window), str(span), outdir, 'unconstraint']
            )  # unconstraint
            outlist.append(
                [
                    fa,
                    fn,
                    gibbs_u,
                    dg_u,
                    nrg_u,
                    printcons,
                    str(window),
                    str(span),
                    outdir,
                    'constraint_unpaired',
                ]
            )  # constraint_paired
            outlist.append(
                [fa, fn, gibbs_p, dg_p, nrg_p, printcons, str(window), str(span), outdir, 'constraint_paired']
            )  # constraint_unpaired

            if os and oe:
                log.debug(
                    logid
                    + 'Second constraint: '
                    + str(','.join(map(str, [goi, data['seq'], len(data['seq']), s, e, os, oe])))
                )
                # enforce both constraints
                # enforce both paired
                fc_p = constrain_paired(fc_p, os, oe)
                # enforce both unpaired
                fc_u = constrain_unpaired(fc_u, os, oe)
                # calculate probs and nrg
                gibbs_u = calc_gibbs(fc_u)
                bppm_u = get_bppm(fc_u.bpp(), cstart, cend)
                dg_u = gibbs_u - gibbs

                gibbs_p = calc_gibbs(fc_p)
                bppm_p = get_bppm(fc_p.bpp(), cstart, cend)
                dg_p = gibbs_p - gibbs

                bpp = calc_bpp(bppm)
                bpp_u = calc_bpp(bppm_u)
                bpp_p = calc_bpp(bppm_p)
                nrg = calc_nrg(bpp)
                nrg_u = calc_nrg(bpp_u)
                nrg_p = calc_nrg(bpp_p)

                outlist.append(
                    [
                        fa,
                        fn,
                        gibbs_u,
                        dg_u,
                        nrg_u,
                        printcons,
                        str(window),
                        str(span),
                        outdir,
                        'bothconstraint_unpaired',
                    ]
                )  # bothconstraint_unpaired
                outlist.append(
                    [
                        fa,
                        fn,
                        gibbs_p,
                        dg_p,
                        nrg_p,
                        printcons,
                        str(window),
                        str(span),
                        outdir,
                        'bothconstraint_paired',
                    ]
                )  # bothconstraint_paired

                # enforce second constraint
                # First clear old constraints
                fc_u = RNA.fold_compound(data['seq'], md)
                fc_p = RNA.fold_compound(data['seq'], md)
                # enforce second paired
                fc_p = constrain_paired(fc_p, os, oe)
                # enforce second unpaired
                fc_u = constrain_unpaired(fc_u, os, oe)
                # calculate probs and nrg
                gibbs_u = calc_gibbs(fc_u)
                bppm_u = get_bppm(fc_u.bpp(), cstart, cend)
                dg_u = gibbs_u - gibbs

                gibbs_p = calc_gibbs(fc_p)
                bppm_p = get_bppm(fc_p.bpp(), cstart, cend)
                dg_p = gibbs_p - gibbs

                bpp = calc_bpp(bppm)
                bpp_u = calc_bpp(bppm_u)
                bpp_p = calc_bpp(bppm_p)
                nrg = calc_nrg(bpp)
                nrg_u = calc_nrg(bpp_u)
                nrg_p = calc_nrg(bpp_p)

                outlist.append(
                    [
                        fa,
                        fn,
                        gibbs_u,
                        dg_u,
                        nrg_u,
                        printcons,
                        str(window),
                        str(span),
                        outdir,
                        'secondconstraint_unpaired',
                    ]
                )  # secondconstraint_unpaired
                outlist.append(
                    [
                        fa,
                        fn,
                        gibbs_p,
                        dg_p,
                        nrg_p,
                        printcons,
                        str(window),
                        str(span),
                        outdir,
                        'secondconstraint_paired',
                    ]
                )  # secondconstraint_paired

        write_out(outlist)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + ''.join(tbe.format()))


def constrain_temp(fa, temp, window, span, an, save, outdir, queue=None, configurer=None, level=None):
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

    logid = SCRIPTNAME + '.constrain_temp: '
    log.info(logid + 'Constraining Temp to ' + temp)
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
        data = {'seq': str(seq_record.seq), 'stru': '', 'nrg': ''}
        # set model details
        md = RNA.md()
        # create new fold_compound object
        fc = RNA.fold_compound(data['seq'], md, RNA.OPTION_PF)
        gibbs = fc.pf()
        bppm = fc.bpp()

        data['stru'] = gibbs[0]
        data['nrg'] = gibbs[1]

        # call prop window calculation
        log.debug(logid + ';'.join(map(str, [gibbs, constrain, conslength])))

        for item in bppm:
            for i in range(int(constrain), int(constrain) + conslength):
                log.debug(bppm.index(item), i, item[i])

        at = up_to_array(data_t['up'], None, len(seq_record.seq))

        diff_nt = an - at

        #####Set temp, and rewrite save for temp
        if save:
            log.info(logid + 'SAVINGTEMP')
            write_temp(fa, str(temp), data_t, diff_nt, str(window), outdir)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + ''.join(tbe.format()))


def foldaround(seq, fc, pos, clength, gibbs, nrg, queue=None, configurer=None, level=None):
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
    Gibbs constraint
    """

    logid = SCRIPTNAME + '.foldaround: '
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
        log.error(logid + ''.join(tbe.format()))


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
    Gibbs unconstraint
    """

    logid = SCRIPTNAME + '.fold_unconstraint: '
    try:
        if queue and level:
            configurer(queue, level)

        # set model details
        md = RNA.md()
        # create new fold_compound object
        fc = RNA.fold_compound(seq, md, RNA.OPTION_PF)
        # call prop window calculation
        gibbs_uc = fc.pf()[1]

        return gibbs_uc
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + ''.join(tbe.format()))


def write_out(resultlist):
    """Print results to file

    Parameters
    ----------
    resultlist: List[str] Contains all results generated by constrain

    Returns
    -------
    Print results to file or STDOUT
    """

    try:
        for result in resultlist:
            fa, fname, gibbs, ddg, nrg, const, window, span, outdir, condition = result
            logid = SCRIPTNAME + '.write_out: '
            goi, chrom, strand = idfromfa(seq_record.id)
            temp_outdir = os.path.join(outdir, goi)

            if fname != 'STDOUT':
                if not os.path.exists(temp_outdir):
                    os.makedirs(temp_outdir)
                if not os.path.exists(
                    os.path.join(
                        temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span]) + '.gz'
                    )
                ):
                    o = gzip.open(
                        os.path.join(
                            temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span]) + '.gz'
                        ),
                        'wb',
                    )
                    o.write(
                        bytes(
                            str.join(
                                '\t', ['Condition', 'FreeNRG(gibbs)', 'deltaG', 'OpeningNRG', 'Constraint']
                            )
                            + '\n',
                            encoding='UTF-8',
                        )
                    )
                    o.write(
                        bytes(
                            str.join('\t', [condition, str(gibbs), str(ddg), str(nrg), str(const)]) + '\n',
                            encoding='UTF-8',
                        )
                    )
                else:
                    log.info(
                        os.path.join(
                            temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span]) + '.gz'
                        )
                        + ' exists, will append!'
                    )
                    o = gzip.open(
                        os.path.join(
                            temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span]) + '.gz'
                        ),
                        'ab',
                    )
                    o.write(
                        bytes(
                            str.join('\t', [condition, str(gibbs), str(ddg), str(nrg), str(const)]) + '\n',
                            encoding='UTF-8',
                        )
                    )
            else:
                print(str.join('\t', ['Condition', 'FreeNRG(gibbs)', 'deltaG', 'OpeningNRG', 'Constraint']))
                print(str.join('\t', [condition, str(gibbs), str(ddg), str(nrg), str(const)]))

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + ''.join(tbe.format()))


def checkexisting(sid, paired, unpaired, cons, window, span, outdir, queue=None, configurer=None, level=None):
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

    logid = SCRIPTNAME + '.checkexisting: '
    try:
        goi, chrom, strand = idfromfa(sid)
        temp_outdir = os.path.join(outdir, goi)

        if os.path.exists(
            os.path.join(temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span]) + '.gz')
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
        log.error(logid + ''.join(tbe.format()))
    return 1


def main():
    """Main process, prepares run_settings dict, creates logging process queue and worker processes for folding, calls fold

    Parameters
    ----------

    Returns
    -------
    Call to fold
    """

    logid = SCRIPTNAME + '.main: '
    try:
        args = parseargs_foldcons()
        queue, listener, worker_configurer = rissmed_logging_setup(args.logdir, args.loglevel, SCRIPTNAME)

        log.info(logid + 'Running ' + SCRIPTNAME + ' on ' + str(args.procs) + ' cores.')
        log.info(
            logid
            + 'CLI: '
            + sys.argv[0]
            + ' '
            + '{}'.format(' '.join([shlex.quote(s) for s in sys.argv[1:]]))
        )

        run_settings, outdir = preprocess(
            args.sequence, args.constrain, args.conslength, args.outdir, args.genes
        )
        fold(
            run_settings,
            outdir,
            args.window,
            args.span,
            args.unconstraint,
            args.unpaired,
            args.paired,
            args.constrain,
            args.conslength,
            args.save,
            args.procs,
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
        log.error(logid + ''.join(tbe.format()))


####################
####    MAIN    ####
####################

if __name__ == '__main__':
    SCRIPTNAME = os.path.basename(__file__).replace('.py', '')
    logid = SCRIPTNAME + '.main: '
    try:
        main()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + ''.join(tbe.format()))

# ConstraintFold.py ends here
