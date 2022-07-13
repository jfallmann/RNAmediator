# FileProcessor.py ---
#
# Filename: FileProcessor.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Fri Aug 21 10:21:43 2020 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Sep  1 10:30:58 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 8
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

import gzip
import inspect

# own
import logging
import os
import re
import sys
import traceback as tb
from collections import defaultdict
from io import StringIO
from typing import DefaultDict, List

####################
# FILE processing  #
####################

try:
    log = logging.getLogger(__name__)  # use module name
    SCRIPTN = os.path.basename(inspect.stack()[-1].filename).replace(".py", "")
    log.debug("LOGGING IN FileProcessor" + str(SCRIPTN) + str(log) + str(log.handlers))
except Exception:
    EXC_TYPE, EXC_VALUE, EXC_TB = sys.exc_info()
    TBE = tb.TracebackException(
        EXC_TYPE,
        EXC_VALUE,
        EXC_TB,
    )
    print("".join(TBE.format()), file=sys.stderr)


def backup(file):
    logid = SCRIPTN + ".backup: "
    try:
        if os.path.exists(file):
            os.rename(file, file + ".bak")
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def parseseq(sequence):
    logid = SCRIPTN + ".parseseq: "
    try:
        if isinstance(sequence, StringIO):
            seq = sequence

        # elif ( isinstance(sequence, str) and sequence == 'random' ):
        #     rand = "\n".join(createrandseq(length, gc, number, alphabet))
        #     seq = StringIO(rand)
        #     o = gzip.open('Random.fa.gz','wb')
        #     o.write(bytes(rand,encoding='UTF-8'))
        #     o.close()

        elif isinstance(sequence, str) and os.path.isfile(sequence):
            if ".gz" in sequence:
                seq = gzip.open(sequence, "rt")
            else:
                seq = open(sequence, "rt")
        else:
            header = ">Seq1:default:nochrom:(.)"
            s = sequence.upper()
            seq = StringIO("{header}\n{s}".format(header=header, s=s))

        return seq

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def idfromfa(fa_id):
    logid = SCRIPTN + ".idfromfa: "
    # goi, chrom, strand = [None, None, None] # not used in the current code
    # fa_id = fa_id.replace("_", "-")
    try:
        goi, chrom = fa_id.split(":")[::2]
        strand = str(fa_id.split(":")[3].split("(")[1][0])
    except (IndexError, ValueError):
        log.warning(
            logid + "Fasta header is not in expected format, you will loose information on strand and chromosome"
        )
        goi = fa_id
        chrom, strand = ["na", "na"]

    if goi and chrom and strand:
        return [str(goi), str(chrom), str(strand)]
    else:
        log.error(logid + "Could not assign any value from fasta header, please check your fasta files")
        sys.exit("Could not assign any value from fasta header, please check your fasta files")


def parse_annotation_bed(bed, annotated=None):
    logid = SCRIPTN + ".parse_annotation_bed: "
    anno = defaultdict(list)
    if os.path.isfile(os.path.abspath(bed)):
        if ".gz" in bed:
            f = gzip.open(os.path.abspath(bed), "rt")
        else:
            f = open(os.path.abspath(bed), "rt")
    else:
        raise FileNotFoundError(f"File {bed} not found")
    try:
        for line in f:
            entries = line.rstrip().split("\t")
            goi = entries[3]
            strand = entries[5]
            if annotated:
                start = int(entries[10]) + 1
                end = int(entries[11])
                strand = entries[14]
            else:
                start = int(entries[1]) + 1
                end = int(entries[2])
            anno[str(goi)].append("|".join(["-".join([str(start), str(end)]), strand]))  # Need strand info here!
        return anno
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def parse_annotation_bed_by_coordinates(bed, annotated=None):
    logid = SCRIPTN + ".parse_annotation_bed_withchrom: "
    anno = defaultdict(list)
    # tmp = defaultdict(list)
    if os.path.isfile(os.path.abspath(bed)):
        if ".gz" in bed:
            f = gzip.open(os.path.abspath(bed), "rt")
        else:
            f = open(os.path.abspath(bed), "rt")
    else:
        raise FileNotFoundError(f"File {bed} not found")
    try:
        for line in f:
            entries = line.rstrip().split("\t")
            chrom = entries[0]
            goi = entries[3]
            strand = entries[5]
            if annotated:
                start = int(entries[10]) + 1
                end = int(entries[11])
                strand = entries[14]
            else:
                start = int(entries[1]) + 1
                end = int(entries[2])
            # tmp[str(goi)].append("|".join(["-".join([str(chrom), str(start), str(end)]), strand]))  # Need strand info here!
            anno[str(goi)].append(
                "|".join(["-".join([str(chrom), str(start), str(end)]), strand])
            )  # Need strand info here! WE ASSUME SORTED, GETS TOO COMPLICATED OTHERWISE
        # We now want to sort by chrom, start to have sorted list in the end, otherwise bigwig will not be generated
        # sortedkeys = sorted(tmp.keys(), key=lambda x: (tmp[x][0].split("|")[0].split("-")[0], tmp[x][0].split("|")[0].split("-")[1]))
        # log.debug(f"{logid} SORTED:{sortedkeys}")
        # for k in sortedkeys:
        #    anno[k] = tmp[k]
        return anno
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def read_constraints_from_bed(bed, linewise=None, ctype="hard"):
    logid = SCRIPTN + ".readConstraintsFromBed: "
    cons = defaultdict(list)
    try:
        for line in bed:
            entries = line.rstrip().split("\t")
            start = int(entries[1]) + 1
            end = entries[2]
            goi = entries[3]
            value = entries[4] if ctype != "hard" else "."
            strand = entries[5]

            if linewise:
                cons["lw"].append("|".join(["-".join([str(start), str(end)]), strand, value]))
            else:
                cons[str(goi)].append("|".join(["-".join([str(start), str(end)]), strand, value]))
        return cons
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def read_paired_constraints_from_bed(bed, linewise=None, constraintype="hard"):
    logid = SCRIPTN + ".readPairedConstraintsFromBed: "
    cons = defaultdict(list)
    try:
        for line in bed:
            entries = line.rstrip().split("\t")
            if len(entries) % 2:
                raise Exception(
                    "Unbalanced paired bed, please make sure the paired bed consists of equal number of "
                    "fields for both constraint entries"
                )
            else:
                second = int((len(entries) / 2) + 1)
            if int(entries[1]) > -1 and int(entries[second]) > -1:
                start_one = int(entries[1]) + 1
                end_one = entries[2]
                goi = entries[3]
                value = entries[4] if constraintype != "hard" else "."
                strand = entries[5]
                start_two = int(entries[second]) + 1
                end_two = int(entries[second + 1])
                value_two = entries[second + 3] if constraintype != "hard" else "."
                if linewise:
                    cons["lw"].append(
                        ":".join(
                            [
                                "|".join(["-".join([str(start_one), str(end_one)]), strand, value]),
                                "|".join(["-".join([str(start_two), str(end_two)]), strand, value_two]),
                            ]
                        )
                    )
                else:
                    cons[str(goi)].append(
                        ":".join(
                            [
                                "|".join(["-".join([str(start_one), str(end_one)]), strand, value]),
                                "|".join(["-".join([str(start_two), str(end_two)]), strand, value_two]),
                            ]
                        )
                    )
        return cons
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def read_constraints_from_csv(csv, linewise=None, constraintype="hard"):
    logid = SCRIPTN + ".readConstraintsCSV: "
    cons: DefaultDict[any, List] = defaultdict(list)
    try:
        for line in csv:
            entries = line.rstrip().split(",")
            start = entries[1]
            end = entries[2]
            goi = entries[3]
            value = entries[4] if constraintype != "hard" else "."
            strand = entries[5]
            if linewise:
                cons["lw"].append("|".join(["-".join([str(start), str(end)]), strand, value]))
            else:
                cons[entries[3]].append("|".join(["-".join([str(start), str(end)]), strand, value]))
        return cons
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def read_constraints_from_generic(generic, linewise=None, constraintype="hard"):
    logid = SCRIPTN + ".readConstraintsFromGeneric: "
    cons: DefaultDict[any, List] = defaultdict(list)

    try:
        for line in generic:
            entries = re.split(r'[ ,|;"]+', line.rstrip())
            if len(entries) > 3:
                goi = entries[1]
                start = entries[2]
                end = entries[3]
                value = entries[4] if constraintype != "hard" else "."
                strand = entries[5]
                if linewise:
                    cons["lw"].append("|".join(["-".join([str(start), str(end)]), strand, value]))
                else:
                    cons[entries[0]].append("|".join(["-".join([str(start), str(end)]), strand]))
            else:
                goi = entries[1]
                start = entries[2]
                end = entries[3]
                value = "."
                strand = "."
                if linewise:
                    cons["lw"].append("|".join(["-".join([str(start), str(end)]), strand, value]))
                else:
                    cons["generic"].append("|".join(["-".join([str(start), str(end)]), strand, value]))
        return cons
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def make_outdir(outdir):
    logid = SCRIPTN + ".makeoutdir: "
    try:
        if not os.path.isabs(outdir):
            outdir = os.path.abspath(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        log.debug(f"{logid} OUTDIR: {outdir}")
        return outdir
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


# TODO: This is in the ConstraintPLFold File isnt it?
# Write results
# def prepare_write_ucons(save, sid, seq, unconstrained, data, region, window, span, outdir, rawentry=None):
#     logid = f"{SCRIPTN}.prepare write_ucons"
#     try:
#         goi, chrom, strand = idfromfa(sid)
#         temp_outdir = os.path.join(outdir, goi)
#     except Exception:
#         exc_type, exc_value, exc_tb = sys.exc_info()
#         tbe = tb.TracebackException(
#             exc_type, exc_value, exc_tb,
#         )
#         log.error(logid + ''.join(tbe.format()))
#
#     try:
#         gr = str(sid.split(':')[3].split('(')[0])
#     except:
#         gr = 'na'
#
#     try:
#         if unconstrained != 'STDOUT':
#             if not os.path.exists(temp_outdir):
#                 os.makedirs(temp_outdir)
#             if rawentry:
#                 if save > 0 and not os.path.exists(os.path.join(temp_outdir, str(
#                         goi + '_' + chrom + '_' + strand + '_' + unconstrained + '_' + rawentry + '_' +
#                         window + '_' + str(
#                             span) + '.gz'))):
#                     with gzip.open(os.path.join(temp_outdir,
#                                                 goi + '_' + chrom + '_' + strand + '_' + unconstrained + '_' +
#                                                 rawentry + '_' + window + '_' + str(
#                                                     span) + '.gz'), 'wb') as o:
#                         out = print_up(data['up'], len(seq), region)
#                         if out and len(out) > 1:
#                             o.write(bytes(out, encoding='UTF-8'))
#                         else:
#                             log.error(logid + "No output produced " + sid)
#             else:
#                 if save > 0 and not os.path.exists(os.path.join(temp_outdir, str(
#                         goi + '_' + chrom + '_' + strand + '_' + unconstrained + '_' + str(
#                             gr) + '_' + window + '_' + str(span) + '.gz'))):
#                     with gzip.open(os.path.join(temp_outdir,
#                                                 goi + '_' + chrom + '_' + strand + '_' + unconstrained + '_' + str(
#                                                     gr) + '_' + window + '_' + str(span) + '.gz'), 'wb') as o:
#                         out = print_up(data['up'], len(seq), region)
#                         if out and len(out) > 1:
#                             o.write(bytes(out, encoding='UTF-8'))
#         else:
#             print(print_up(data, len(seq), region))
#     except Exception:
#         exc_type, exc_value, exc_tb = sys.exc_info()
#         tbe = tb.TracebackException(
#             exc_type, exc_value, exc_tb,
#         )
#         log.error(logid + ''.join(tbe.format()))
#     return 1
#
#
# def prepare_write_cons(save, sid, seq, paired, unpaired, data_u, data_p, constrain, region, diff_nu, diff_np, window,
#                        span, outdir):
#     try:
#         goi, chrom, strand = idfromfa(sid)
#         temp_outdir = os.path.join(outdir, goi)
#     except Exception:
#         exc_type, exc_value, exc_tb = sys.exc_info()
#         tbe = tb.TracebackException(
#             exc_type, exc_value, exc_tb,
#         )
#         log.error(logid + ''.join(tbe.format()))
#     # print outputs to file or STDERR
#     try:
#         if paired != 'STDOUT':
#             if not os.path.exists(temp_outdir):
#                 os.makedirs(temp_outdir)
#             if save > 0 and not os.path.exists(os.path.join(temp_outdir,
#                                                             'StruCons_' + goi + '_' + chrom + '_' + strand + '_' +
#                                                             constrain + '_' + paired + '_' + window + '_' + str(
#                                                                 span) + '.gz')):
#                 with gzip.open(os.path.join(temp_outdir,
#                                             'StruCons_' + goi + '_' + chrom + '_' + strand + '_' + constrain + '_' +
#                                             paired + '_' + window + '_' + str(
#                                                 span) + '.gz'), 'wb') as o:
#                     out = print_up(data_p['up'], len(seq), region)
#                     if out and len(out) > 1:
#                         o.write(bytes(out, encoding='UTF-8'))
#                     else:
#                         log.error(logid + "No output produced " + sid)
#         else:
#             print(print_up(data_p['up'], len(seq), region))
#
#         if unpaired != 'STDOUT':
#             if not os.path.exists(temp_outdir):
#                 os.makedirs(temp_outdir)
#             if save > 0 and not os.path.exists(os.path.join(temp_outdir,
#                                                             'StruCons_' + goi + '_' + chrom + '_' + strand + '_' +
#                                                             constrain + '_' + unpaired + '_' + window + '_' + str(
#                                                                 span) + '.gz')):
#                 with gzip.open(os.path.join(temp_outdir,
#                                             'StruCons_' + goi + '_' + chrom + '_' + strand + '_' +
#                                             constrain + '_' + unpaired + '_' + window + '_' + str(
#                                                 span) + '.gz'), 'wb') as o:
#                     out = print_up(data_u['up'], len(seq), region)
#                     if out and len(out) > 1:
#                         o.write(bytes(out, encoding='UTF-8'))
#                     else:
#                         log.warning("No output produced " + sid)
#         else:
#             print(print_up(data_u['up'], len(seq), region))
#
#         if diff_nu.any():
#             if unpaired != 'STDOUT':
#                 if not os.path.exists(temp_outdir):
#                     os.makedirs(temp_outdir)
#                 if not os.path.exists(os.path.join(temp_outdir,
#                                                    'StruCons_' + goi + '_' + chrom + '_' + strand + '_' + constrain +
#                                                    '_diffnu_' + window + '_' + str(
#                                                        span) + '.npy')):
#                     printdiff(diff_nu, os.path.join(temp_outdir,
#                                                     'StruCons_' + goi + '_' + chrom + '_' + strand + '_' + constrain +
#                                                     '_diffnu_' + window + '_' + str(
#                                                         span) + '.npy'))
#             else:
#                 npprint(diff_nu)
#
#         if diff_np.any():
#             if unpaired != 'STDOUT':
#                 if not os.path.exists(temp_outdir):
#                     os.makedirs(temp_outdir)
#                 if not os.path.exists(os.path.join(temp_outdir,
#                                                    'StruCons_' + goi + '_' + chrom + '_' + strand + '_' + constrain +
#                                                    '_diffnp_' + window + '_' + str(
#                                                        span) + '.npy')):
#                     printdiff(diff_np, os.path.join(temp_outdir,
#                                                     'StruCons_' + goi + '_' + chrom + '_' + strand + '_' + constrain +
#                                                     '_diffnp_' + window + '_' + str(
#                                                         span) + '.npy'))
#             else:
#                 npprint(diff_np)
#     except Exception:
#         exc_type, exc_value, exc_tb = sys.exc_info()
#         tbe = tb.TracebackException(
#             exc_type, exc_value, exc_tb,
#         )
#         log.error(logid + ''.join(tbe.format()))
#     return 1

#
# FileProcessor.py ends here
