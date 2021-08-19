from __future__ import annotations  # It will become the default in Python 3.10.
from typing import Iterable, Tuple, Union, Dict, List
from dataclasses import dataclass
import logging
import datetime
import os
import sys
import traceback as tb
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from RIssmed.Tweaks.FileProcessor import (
    parseseq,
    idfromfa,
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
import multiprocessing


log = logging.getLogger(__name__)  # use module name
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
     constrainlist: Iterable[Tuple[Constraint]], optional
         List of constraint objects. (default is None)
     genomic_coords: Constraint, optional optional
         genomic coordinates of the sequence record (default is None)
    """

    def __init__(
        self,
        sequence_record: SeqIO.SeqRecord,
        gene: str = "nogene",
        chrom: str = "nochrom",
        strand: str = "+",
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

        self._check_strands()

    def _check_strands(self):
        if self._constrainlist is not None:
            for constraint_tuple in self._constrainlist:
                for entry in constraint_tuple:
                    if self.strand in ["+", "-"]:
                        assert (
                            entry.strand == self.strand
                        ), "strand values of constraint does not match the strand from the sequence"
                    assert (
                        entry.__class__ == Constraint
                    ), "can only add Contraint objects to the constraintlist"

    def add_constraints(self, constraints: Tuple[Constraint]):
        """adds a constraints to the list of constraints"""
        for constraint in constraints:
            assert (
                constraint.__class__ == Constraint
            ), "can only add Contraint objects to the constraintliste"
            if self.strand in ["+", "-"]:
                assert (
                    self.strand == constraint.strand
                ), "strand values of constraint does not match the strand from the sequence"
        self._constrainlist += [constraints]

    @property
    def constrainlist(self):
        """getter method for constraintlist attribute"""
        return self._constrainlist


@dataclass(frozen=True)
class Constraint:
    start: int
    end: int
    strand: str

    def __str__(self):
        return f"{self.start}-{self.end}|{self.strand}"


def get_gene_coords(
    genecoords: Union[None, Dict], goi: str, strand: str
) -> Tuple[int, int, str]:
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
    Tuple[int, int, str]
        genomic start, end, strand retrieved from genecoords dict or a default value (0, 0, '.')
    """
    logid = f"{SCRIPTNAME}.read_constraints "
    if genecoords:
        if goi in genecoords:
            gs, ge, gstrand = get_location(genecoords[goi][0])
            if gstrand != strand:
                log.warning(
                    logid
                    + "Strand values differ between Gene annotation and FASTA file! Please check your input for "
                    + str(goi)
                )
        else:
            gs, ge, gstrand = 0, 0, "."
            log.warning(
                logid
                + "No coords found for gene "
                + goi
                + "! Assuming coordinates are already local!"
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
    return gs, ge, gstrand


def set_run_settings_dict(
    sequence, constrain: str, conslength: int, genes: str
) -> Dict[str, SequenceSettings]:
    """Use command line parameters to build the run settings dictionary.

    Parameters
    ----------
     sequence : str
        The file location of the sequence
     constrain : str
        The file location of constrain file
     conslength : int
         Length of the constraint, only used if constrain is sliding
     genes:
         The file location of the genomic coordinates bed file

    Returns
    -------
    Dict[str, SequenceSettings]
        a dictionary using the fasta sequence id as key and stores corresponding settings in an SequenceSettings
        object
    """
    run_settings: Dict[str, SequenceSettings] = dict()
    sequence = parseseq(sequence)
    if "ono" == str(constrain.split(",")[0]):
        constrain = constrain.split(",")[1]
        constraintlist = read_constraints(constrain, linewise=True)
        for x, record in enumerate(SeqIO.parse(sequence, "fasta")):
            goi, chrom, strand = idfromfa(record.id)
            cons = constraintlist["lw"][x]
            run_settings = add_rissmed_constraint(
                run_settings, cons, record, goi, chrom, strand
            )
    elif constrain == "sliding":
        for record in SeqIO.parse(sequence, "fasta"):
            goi, chrom, strand = idfromfa(record.id)
            for start in range(1, len(record.seq) - conslength + 2):
                end = start + conslength - 1
                cons = str(start) + "-" + str(end) + "|" + str(strand)  # AUA
                run_settings = add_rissmed_constraint(
                    run_settings, cons, record, goi, chrom, strand
                )
    else:
        constraintlist = read_constraints(constrain=constrain)
        for x, record in enumerate(SeqIO.parse(sequence, "fasta")):
            goi, chrom, strand = idfromfa(record.id)
            cons = (
                constraintlist[goi] if type(constraintlist) == defaultdict else constraintlist
            )
            for entry in cons:
                run_settings = add_rissmed_constraint(
                    run_settings, entry, record, goi, chrom, strand
                )
    if genes != "":
        # get genomic coords to print to bed later, should always be just one set of coords per gene
        genecoords = parse_annotation_bed(genes)
    else:
        genecoords = None
    for entry in run_settings:
        fasta_settings = run_settings[entry]
        goi, chrom, strand = idfromfa(fasta_settings.sequence_record.id)
        genomic_start, genomic_end, genomic_strand = get_gene_coords(
            genecoords, goi, strand
        )
        fasta_settings.genomic_coords = Constraint(
            genomic_start, genomic_end, genomic_strand
        )

    return run_settings


def add_rissmed_constraint(
    run_settings: Dict[str, SequenceSettings],
    constraints: str,
    record: SeqIO.SeqRecord,
    goi: str = "nogene",
    chrom: str = "nochrom",
    sequence_strand: str = "+",
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

    Returns
    -------
    Dict[str, SequenceSettings]
        a dictionary using the fasta sequence id as key and stores corresponding settings in an SequenceSettings
        object.
    """

    cons_list = []
    for constraint in constraints.split(
        ":"
    ):  # Should now work with paired constraints split via : separator
        cons = constraint.split("|")
        cons_strand = cons[1] if len(cons) > 1 else sequence_strand
        cons = cons[0]
        cons_start, cons_end = cons.split("-")
        cons_list.append(Constraint(int(cons_start), int(cons_end), cons_strand))
    cons_tuple = tuple(cons_list)
    if record.id in run_settings:
        run_settings[record.id].add_constraints(cons_tuple)
    else:
        settings = SequenceSettings(
            record,
            constrainlist=[cons_tuple],
            chrom=chrom,
            gene=goi,
            strand=sequence_strand,
        )
        run_settings[record.id] = settings
    return run_settings


# put this into Fileprocessing ?
def read_constraints(constrain: str, linewise: bool = False) -> Dict[str, List[str]]:
    """Reads constraints from the constraints file

    Parameters
    ----------
     constrain : Union[None, Dict]
        file location of the constrains file
     linewise : bool, optional
        does not add the gene identifier to the returned cinstrantslist dictionary is set to True
        (default is False)

    Returns
    -------
    Dict[str, List[str]]
        Dictionary containing the goi in the bed file as keys and constraints in string format (start-end|strand)
        as values. If linewise is set to True only a single key ("lw") is used to store all constraints.
    """
    constraintlist = []
    logid = f"{SCRIPTNAME}.read_constraints"
    if os.path.isfile(constrain):
        if ".bed" in constrain:
            log.info(logid + "Parsing constraints from Bed " + constrain)
            if ".gz" in constrain:
                f = gzip.open(constrain, "rt")
            else:
                f = open(constrain, "rt")
            if "paired" in constrain:
                constraintlist = read_paired_constraints_from_bed(
                    f, linewise
                )  # Not sure if it works but it should
            else:
                constraintlist = read_constraints_from_bed(f, linewise)
        elif ".csv" in constrain:
            if ".gz" in constrain:
                f = gzip.open(constrain, "rt")
            else:
                f = open(constrain, "rt")
            constraintlist = read_constraints_from_csv(f, linewise)
        else:
            if ".gz" in constrain:
                f = gzip.open(constrain, "rt")
            else:
                f = open(constrain, "rt")
            constraintlist = read_constraints_from_generic(f, linewise)
        f.close()
    elif constrain == "file" or constrain == "paired":
        log.info(
            logid
            + "Calculating probs for constraint from file "
            + str(goi + "_constraints")
        )
        with open(goi + "_constraints", "rt") as o:
            for line in o:
                conslist.append(line.rstrip())
    elif constrain == "none":
        constraintlist = ["NOCONS"]
    elif constrain == "sliding":
        constraintlist = list()
    elif "-" in constrain:
        log.info(logid + "Calculating probs for constraint " + constrain)
        if linewise is False:
            constraintlist = constrain.split(",")
        else:
            constraintlist = defaultdict(list)
            for cons in constrain.split(","):
                constraintlist["lw"].append(cons)

    elif constrain == "temperature":
        log.info(logid + "Calculating probs for temperature constraint" + temprange)
        raise NotImplementedError("Temperature range folding needs to be reimplemented")
    else:
        log.error(logid + "Could not compute constraints from input " + str(constrain))
        sys.exit()
    return constraintlist


def preprocess(sequence: str, constrain: str, conslength: int, outdir: str, genes: str):
    """builds the run settings dict and creates the output directory

    Parameters
    ----------
     sequence : str
        The file location of the sequence
     constrain : str
        The file location of constrain file
     conslength : int
         Length of the constraint, only used if constrain is sliding
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
            if not os.path.isabs(outdir):
                outdir = os.path.abspath(outdir)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
        else:
            outdir = os.path.abspath(os.getcwd())

        run_settings = set_run_settings_dict(sequence, constrain, conslength, genes)

        return run_settings, outdir

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


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
    listener = multiprocessing.Process(
        target=listener_process, args=(queue, listener_configurer, logfile, loglevel)
    )
    listener.start()
    worker_configurer(queue, loglevel)

    return queue, listener, worker_configurer


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


def expand_window(start, end, window, multiplyer, seqlen):
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


def localize_window(start, end, window, seqlen, multiplyer=2):
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
