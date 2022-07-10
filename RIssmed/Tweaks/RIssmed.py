from __future__ import annotations
from typing import Iterable, Tuple, Union, Dict, List
from dataclasses import dataclass
import logging
import datetime
import os
import sys
import traceback as tb
from numpy.random import default_rng
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from RIssmed.Tweaks.FileProcessor import (
    parseseq,
    idfromfa,
    make_outdir,
    parse_annotation_bed,
    read_constraints_from_bed,
    read_constraints_from_csv,
    read_constraints_from_generic,
    read_paired_constraints_from_bed,
)
import gzip
from RIssmed.Tweaks.RNAtweaks import get_location
from RIssmed.Tweaks.logger import (
    makelogdir,
    makelogfile,
    listener_process,
    listener_configurer,
    worker_configurer,
)
from RIssmed.Tweaks.Collection import check_run
import multiprocessing


log = logging.getLogger(__name__)  # use module name
# log.propagate = True
SCRIPTNAME = os.path.basename(__file__).replace(".py", "")


class SequenceSettings:
    """Constraint(PL)fold settings for a sequence object.

     Attributes
    ----------
     sequence_record : SeqIO.SeqRecord
        Sequence record the settings belong to
     gene : str, optional
        name of the gene of interest (default is nogene)
     chrom: str, optional
         chromosome of the gene (default is nochrom)
     strand: str, optional
         strand of the gene (default is +)
    constraintype: str, optional
        Type of constraint to apply, can be ['hard'(Default), 'soft'(not implemented yet), 'mutate']
     constrainlist: Iterable[Tuple[Constraint]], optional
         List of constraint objects. (default is None)
     genomic_coords: Constraint, optional optional
         genomic coordinates of the sequence record (default is None)

    """

    @check_run
    def __init__(
        self,
        sequence_record: SeqIO.SeqRecord,
        gene: str = "nogene",
        chrom: str = "nochrom",
        strand: str = "+",
        constraintype: str = "hard",
        constrainlist: Iterable[Tuple[Constraint]] = None,
        genomic_coords: Constraint = None,
    ):
        sequence_record.seq = Seq(
            str(sequence_record.seq).upper().replace("T", "U")
        )  # We always want RNA Sequence to have consistent output to ViennaRNA-CLI
        self.sequence_record = sequence_record
        self._constrainlist = list(constrainlist) if constrainlist is not None else []
        if strand in ["+", "-", "na", "."]:
            self.strand = strand
        else:
            self.strand = "na"
            log.warning("strand value automatically set to na")
        self.chromosome = chrom
        self.gene = gene
        self.genomic_coords = genomic_coords
        self._constraintype = constraintype

        self._check_strands()

    @check_run
    def _check_strands(self):
        if self._constrainlist is not None:
            for constraint_tuple in self._constrainlist:
                for entry in constraint_tuple:
                    if entry is not None:
                        if self.strand in ["+", "-"]:
                            assert (
                                entry.strand == self.strand
                            ), "strand values of constraint does not match the strand from the sequence"
                        assert entry.__class__ == Constraint, "can only add Contraint objects to the constraintlist"

    @check_run
    def add_constraints(self, constraints: Tuple[Constraint]):
        """adds a constraint to the list of constraints"""
        for constraint in constraints:
            assert constraint.__class__ == Constraint, "can only add Contraint objects to the constraintliste"
            if self.strand in ["+", "-"]:
                try:
                    assert self.strand == constraint.strand
                except AssertionError:
                    print("AssertionError: strand values of constraint does not match the strand from the sequence")
                    raise

        self._constrainlist += [constraints]

    @check_run
    def set_constraintype(self, ctype: str):
        """sets constraint type"""
        self._constraintype = ctype

    @property
    def constrainlist(self):
        """getter method for constraintlist attribute"""
        return self._constrainlist

    @property
    def constraintype(self):
        """getter method for constraintype attribute"""
        return self._constraintype


@dataclass(frozen=True)
class Constraint:
    start: int
    end: int
    strand: str
    type: str = ""

    def __str__(self):
        if self.type != "":
            return f"{self.start}-{self.end}|{self.strand}|{self.type}"
        else:
            return f"{self.start}-{self.end}|{self.strand}"


@check_run
def get_gene_coords(genecoords: Union[None, Dict], goi: str, strand: str) -> Tuple[int, int, str, str]:
    """Get genomic coordinates for a gene (goi) from the genecoords dict or returns  default values

    Parameters
    ----------
     genecoords : Union[None, Dict]
        gene coords dictionary, storing start, end and strand of a gene (key) as values
     goi : str
        gene of interest, for which coordinates should be retrieved
     strand:
         strand of the gene to compare whether strand matches the strand in the gencoords dict

    Returns
    -------
    Tuple[int, int, str, str]
        genomic start, end, strand retrieved from genecoords dict or a default value (0, 0, '.')
    """
    logid = f"{SCRIPTNAME}.get_gene_coords "
    log.debug(f"{logid} Goi:{goi}, coords:{genecoords}")
    if genecoords:
        if goi in genecoords:
            gs, ge, gstrand, value = get_location(genecoords[goi][0])
            if gstrand != strand:
                log.warning(
                    logid
                    + "Strand values differ between Gene annotation and FASTA file! Please check your input for "
                    + str(goi)
                )
        else:
            gs, ge, gstrand, value = 0, 0, ".", "."
            log.warning(logid + "No coords found for gene " + goi + "! Assuming coordinates are already local!")
    else:
        gs = ge = 0
        gstrand = "."
        value = "."
        log.warning(logid + "No coords found for gene " + goi + "! Assuming coordinates are already local!")
    return gs, ge, gstrand, value


@check_run
def set_run_settings_dict(
    sequence, constraint: str, conslength: int, genes: str, constraintype: str = "hard"
) -> Dict[str, SequenceSettings]:
    """Use command line parameters to build the run settings dictionary.

    Parameters
    ----------
     sequence : str
        The file location of the sequence
     constraint : str
        The file location of constrain file
     conslength : int
        Length of the constraint, only used if constrain is sliding
    genes:
        The file location of the genomic coordinates bed file
     constraintype : str, optional
        Type of constraint to apply, can be ['hard'(Default), 'soft'(not implemented yet), 'mutate']

    Returns
    -------
    Dict[str, SequenceSettings]
        a dictionary using the fasta sequence id as key and stores corresponding settings in an SequenceSettings
        object
    """
    logid = SCRIPTNAME + ".set_run_settings_dict: "
    log.debug(f"{logid} seq:{sequence} constraint:{constraint} genes:{genes} constraintype:{constraintype}")
    try:
        run_settings: Dict[str, SequenceSettings] = dict()
        sequence = parseseq(sequence)
        if genes != "":
            # get genomic coords to print to bed later, should always be just one set of coords per gene
            genecoords = parse_annotation_bed(genes)
        else:
            genecoords = None
        if "ono" == str(constraint.split(",")[0]):
            constraint = constraint.split(",")[1]  # 0-based
            constraintlist = read_constraints(constraint, linewise=True, constraintype=constraintype)
            for x, record in enumerate(SeqIO.parse(sequence, "fasta")):
                goi, chrom, strand = idfromfa(record.id)
                cons = constraintlist["lw"][x]
                run_settings = add_rissmed_constraint(
                    run_settings, cons, record, goi, chrom, strand, constraintype=constraintype
                )
        elif constraint == "sliding":
            for record in SeqIO.parse(sequence, "fasta"):
                goi, chrom, strand = idfromfa(record.id)
                for start in range(1, len(record.seq) - conslength + 2):  # 0-based
                    end = start + conslength - 1
                    cons = str(start) + "-" + str(end) + "|" + str(strand)
                    run_settings = add_rissmed_constraint(
                        run_settings, cons, record, goi, chrom, strand, constraintype=constraintype
                    )
        elif constraint[:7] == "random,":
            if constraintype != "hard":
                raise NotImplementedError("Random constraints are currently only supported for hard constraints")

            nr_cons = int(constraint.split(",")[1]) - 1
            rng = default_rng()
            for record in SeqIO.parse(sequence, "fasta"):
                goi, chrom, strand = idfromfa(record.id)
                genomic_start, genomic_end, genomic_strand, value = get_gene_coords(genecoords, goi, strand)
                log.debug(f"{logid} CHECKCOORDS: {genomic_start}, {genomic_end}, {genomic_strand}, {value}")
                randstarts = rng.integers(low=genomic_start, high=genomic_end - conslength + 2, size=nr_cons)
                for start in randstarts:  # 0-based
                    end = start + conslength - 1
                    cons = str(start) + "-" + str(end) + "|" + str(strand)
                    run_settings = add_rissmed_constraint(
                        run_settings, cons, record, goi, chrom, strand, constraintype=constraintype
                    )
        else:
            constraintlist = read_constraints(constraint=constraint, constraintype=constraintype)
            for x, record in enumerate(SeqIO.parse(sequence, "fasta")):
                goi, chrom, strand = idfromfa(record.id)
                cons = constraintlist[goi] if type(constraintlist) == defaultdict else constraintlist
                log.debug(f"{logid} Setting {goi} {chrom} {strand} constraint {cons}")
                for entry in cons:
                    run_settings = add_rissmed_constraint(
                        run_settings, entry, record, goi, chrom, strand, constraintype=constraintype
                    )
        for entry in run_settings:
            fasta_settings = run_settings[entry]
            goi, chrom, strand = idfromfa(fasta_settings.sequence_record.id)
            genomic_start, genomic_end, genomic_strand, value = get_gene_coords(genecoords, goi, strand)
            fasta_settings.genomic_coords = Constraint(genomic_start, genomic_end, genomic_strand, value)

        return run_settings

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


@check_run
def add_rissmed_constraint(
    run_settings: Dict[str, SequenceSettings],
    constraints: str,
    record: SeqIO.SeqRecord,
    goi: str = "nogene",
    chrom: str = "nochrom",
    sequence_strand: str = "+",
    constraintype: str = "hard",
):
    """Adds constraints in string format to the goi in the run settings dict. Creates sequence settings if missing

    Parameters
    ----------
     run_settings : Dict[str, SequenceSettings]
        run settings dictionary using fasta ids as keys and Sequence Settings as values
     constraints : str
        constraint in string format start-end|strand
     record : SeqIO.SeqRecord
         Sequence record of the sequence to which the constraint should be added
     goi:
         gene name assumed from fasta sequence header
     chrom:
         chromosome of the sequence
     sequence_strand:
        strand of the sequence
     constraintype: str, optional
        Type of constraint to apply, can be ['hard'(Default), 'soft'(not implemented yet), 'mutate']

    Returns
    -------
    Dict[str, SequenceSettings]
        a dictionary using the fasta sequence id as key and stores corresponding settings in an SequenceSettings
        object.
    """

    logid = f"{SCRIPTNAME}.add_rissmed_constraint"
    cons_list = []
    for constraint in constraints.split(":"):  # Should now work with paired constraints split via : separator
        cons = constraint.split("|")
        log.debug(f"{logid} Constraint: {cons}")
        if cons[0] == "NOCONS":
            # cons_list.append(NoConstraint("NOCONS", "NOCONS", sequence_strand))
            cons_list.append(None)
        else:
            cons_strand = cons[1] if len(cons) > 1 else sequence_strand
            cons_type_value = cons[2] if len(cons) > 2 else "hard"
            cons = cons[0]
            cons_start, cons_end = cons.split("-")
            cons_list.append(Constraint(int(cons_start), int(cons_end), cons_strand, cons_type_value))
    cons_tuple = tuple(cons_list)
    if record.id in run_settings:
        run_settings[record.id].add_constraints(cons_tuple)
    else:
        settings = SequenceSettings(
            record,
            constrainlist=[cons_tuple],
            constraintype=constraintype,
            chrom=chrom,
            gene=goi,
            strand=sequence_strand,
        )
        run_settings[record.id] = settings
    return run_settings


# put this into Fileprocessing ?
@check_run
def read_constraints(constraint: str, linewise: bool = False, constraintype: str = "hard") -> Dict[str, List[str]]:
    """Reads constraints from the constraints file

    Parameters
    ----------
    constraint : Union[None, Dict]
        file location of the constrains file
    linewise : bool, optional
        does not add the gene identifier to the returned cinstrantslist dictionary is set to True
        (default is False)
    constraintype : str, optional
        sets type of constraint (default is hard)

    Returns
    -------
    Dict[str, List[str]]
        Dictionary containing the goi in the bed file as keys and constraints in string format (start-end|strand)
        as values. If linewise is set to True only a single key ("lw") is used to store all constraints.
    """
    constraintlist = []
    logid = f"{SCRIPTNAME}.read_constraints: "
    log.debug(f"{logid} constraint:{constraint}, linewise:{linewise} constype:{constraintype}")
    if os.path.isfile(constraint):
        if ".bed" in constraint:
            log.info(logid + "Parsing constraints from Bed " + constraint)
            if ".gz" in constraint:
                f = gzip.open(constraint, "rt")
            else:
                f = open(constraint, "rt")
            if "paired" in constraint:
                constraintlist = read_paired_constraints_from_bed(f, linewise, constraintype)
            else:
                constraintlist = read_constraints_from_bed(f, linewise, constraintype)
        elif ".csv" in constraint:
            if ".gz" in constraint:
                f = gzip.open(constraint, "rt")
            else:
                f = open(constraint, "rt")
            constraintlist = read_constraints_from_csv(f, linewise, constraintype)
        else:
            if ".gz" in constraint:
                f = gzip.open(constraint, "rt")
            else:
                f = open(constraint, "rt")
            constraintlist = read_constraints_from_generic(f, linewise, constraintype)
        f.close()
    # elif constrain == "file" or constrain == "paired":
    #    log.info(
    #        logid
    #        + "Calculating probs for constraint from file "
    #        + str(goi + "_constraints")
    #    )
    #    with open(goi + "_constraints", "rt") as o:
    #        for line in o:
    #            conslist.append(line.rstrip())
    elif constraint == "scanning":
        constraintlist = ["NOCONS"]
    elif constraint == "sliding":
        constraintlist = list()
    elif "-" in constraint:
        log.info(logid + "Calculating probs for constraint " + constraint)
        if linewise is False:
            constraintlist = list()
            for x in constraint.split(","):
                a = x.split("|")
                if len(a) < 2:
                    a.extend([".", "."])
                elif len(a) < 3:
                    if a[1] in ["+", "-", "."]:
                        a.append(".")
                    elif type(a[1]) == str:
                        a.append(a[1])
                        a[1] = "."
                    else:
                        a[1] = a[2] = "."
                x = "|".join(a)
                constraintlist.append(x)
        else:
            constraintlist = defaultdict(list)
            for x in constraint.split(","):
                a = x.split("|")
                for i in [2, 3]:
                    if len(a) < 2:
                        a.extend([".", "."])
                    elif len(a) < 3:
                        if a[1] in ["+", "-", "."]:
                            a.append(".")
                        elif type(a[1]) == str:
                            a.append(a[1])
                            a[1] = "."
                        else:
                            a[1] = a[2] = "."
                x = "|".join(a)
                constraintlist["lw"].append(x)

    elif constraint == "temperature":
        raise NotImplementedError("Temperature range folding needs to be reimplemented")
    else:
        log.error(logid + "Could not compute constraints from input " + str(constraint))
        sys.exit()
    log.debug(f"{logid} Constraintlist: {constraintlist}")
    return constraintlist


@check_run
def preprocess(sequence: str, constraint: str, conslength: int, constype: str, outdir: str, genes: str):
    """builds the run settings dict and creates the output directory

    Parameters
    ----------
     sequence : str
        The file location of the sequence
     constraint : str
        The file location of constrain file
     conslength : int
         Length of the constraint, only used if constrain is sliding
     constype : str
         Type of constraint to apply, can be ['hard'(Default), 'soft'(not implemented yet), 'mutate']
     outdir : str
         Location of the Outpu directory. If it is an empty string os.cwd() is used
     genes:
         The file location of the genomic coordinates bed file

    Returns
    -------
    Tuple[Dict[str, SequenceSettings], str]
        a dictionary using the fasta sequence id as key and stores corresponding settings in an SequenceSettings
        object and the absolute path to the output directory

    """
    logid = SCRIPTNAME + ".preprocess: "
    try:
        # set path for output
        if outdir:
            outdir = make_outdir(outdir)
        else:
            outdir = os.path.abspath(os.getcwd())

        run_settings = set_run_settings_dict(sequence, constraint, conslength, genes, constype)

        return run_settings, outdir

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


@check_run
def rissmed_logging_setup(logdir: str, loglevel: str, runscript: str):
    """creates log process, queue, config, dir and file according to runscript

    Parameters
    ----------
     logdir : str
        The directory for log files
     loglevel : str
        The level for logging
     runscript : str
         Name of the script for which logging is set up

    Returns
    -------
    Tuple[queue, listener, worker_configurer]
        multiprocessing queue, listener process and settings for worker thread logging
    """

    ts = str(datetime.datetime.now().strftime("%Y%m%d_%H_%M_%S_%f"))
    logfile = str.join(os.sep, [os.path.abspath(logdir), runscript + "_" + ts + ".log"])
    makelogdir(logdir)
    makelogfile(logfile)

    if multiprocessing.get_start_method(allow_none=True) != "spawn":
        multiprocessing.set_start_method("spawn")
    queue = multiprocessing.Manager().Queue(-1)
    listener = multiprocessing.Process(target=listener_process, args=(queue, listener_configurer, logfile, loglevel))
    listener.start()
    worker_configurer(queue, loglevel)

    return queue, listener, worker_configurer


@check_run
def expand_pl_window(start, end, window, multiplyer, seqlen):
    logid = SCRIPTNAME + ".expand_window: "
    try:
        tostart = start - multiplyer * window
        if tostart < 1:
            tostart = 1
        toend = end + multiplyer * window
        if toend > seqlen:
            toend = seqlen
        return [tostart, toend]
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


@check_run
def localize_pl_window(start, end, window, seqlen, multiplyer=2):
    logid = SCRIPTNAME + ".localize_window: "
    try:
        diff = start - window
        if diff < 1:
            locws = 1
        else:
            locws = diff
        # this makes sure that if the start was trimmed, we do not just extend too much
        locwe = diff + multiplyer * window + (end - start)

        if locwe > seqlen:
            locwe = seqlen
        return [locws, locwe]
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


@check_run
def expand_window(start, end, window, seqlen=None, dist=None):
    logid = SCRIPTNAME + ".expand_window: "
    try:
        if dist:
            tostart = start - (window - abs(dist))
            toend = end + (window - abs(dist)) + 1
        else:
            tostart = start - window
            toend = end + window + 1

        if tostart < 0:
            tostart = 1
        if toend > seqlen:
            toend = seqlen
        return [tostart, toend]
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


@check_run
def localize_window(start, end, tostart, toend, seqlen=None):
    logid = SCRIPTNAME + ".localize_window: "
    try:
        wstart = start - tostart
        wend = end - tostart
        if wstart < 0 or wend > seqlen:
            log.error(
                logid
                + "start of constraint "
                + str(wstart)
                + " end of constraint "
                + str(wend)
                + " while length of sequence "
                + str(seqlen)
            )
        if wstart < 0:
            wstart = 0
        if wend > seqlen:
            wend = seqlen
        return [wstart, wend]
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))
