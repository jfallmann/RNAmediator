### RIssmedArgparsers.py ---
##
## Filename: RIssmedArgparsers.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Sep  6 09:02:18 2018 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Wed Dec 16 13:19:45 2020 (+0100)
##           By: Joerg Fallmann
##     Update #: 384
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
import os
import sys
import logging

## other modules
import traceback as tb
import argparse
import numpy as np

############################################################
######################## FUNCTIONS #########################
############################################################

try:
    log = logging.getLogger(__name__)  # use module name
    scriptn = __name__  # os.path.basename(inspect.stack()[-1].filename).replace('.py', '')
    log.debug("LOGGING IN Argparsers" + str(scriptn) + str(log) + str(log.handlers))
except Exception:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type,
        exc_value,
        exc_tb,
    )
    print("".join(tbe.format()), file=sys.stderr)


##############################
########## ARGPARSE ##########
##############################


def parseargs_plcons():
    parser = argparse.ArgumentParser(
        description="Calculate base pairing probabilities with and without constraints for sequences with given window size, span and region."
    )
    parser.add_argument("-s", "--sequence", type=str, help="Sequence to fold")
    parser.add_argument("-w", "--window", type=int, default=240, help="Size of window")
    parser.add_argument("-l", "--span", type=int, default=60, help="Length of bp span")
    parser.add_argument("-u", "--region", type=int, default=1, help="Length of region")
    parser.add_argument("-m", "--multi", type=int, default=2, help="Multiplyer for window expansion")
    parser.add_argument(
        "-r",
        "--unconstraint",
        type=str,
        default="raw",
        help="Print output of unconstraint folding to file with this name",
    )
    parser.add_argument(
        "-n",
        "--unpaired",
        type=str,
        default="unpaired",
        help="Print output of unpaired folding to file with this name",
    )
    parser.add_argument(
        "-p",
        "--paired",
        type=str,
        default="paired",
        help="Print output of paired folding to file with this name",
    )
    parser.add_argument(
        "-x",
        "--constrain",
        type=str,
        default="sliding",
        help="Region to constrain, either sliding window (default) or region to constrain (e.g. 1-10) or path to file containing constraint regions, can be BED file or the string 'file'. Latter requires that a file following the naming pattern $fastaID_constraints e.g. Sequence1_constraints can be found. If 'paired', the first entry of the constraint file will become a fixed constraint and paired with all the others, If 'ono' one line of constraint file is used per line of sequence in the sequence file. Choose 'scanning' to fold sequence without constraint. choices = [scanning, sliding, the string 'file' or a filename, paired, or simply 1-10,2-11 or 1-10;15-20,2-11;16-21 for paired or ono(oneonone),filename to use one line of constraint file for one sequence from fasta]",
    )
    parser.add_argument(
        "-y",
        "--conslength",
        type=int,
        default=1,
        help="Length of region to constrain for sliding window",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        type=int,
        default=37,
        help="Temperature for structure prediction",
    )
    parser.add_argument(
        "--save",
        type=int,
        default=0,
        choices=[0, 1],
        help="Save the output as numpy files only [0] or also as gzipped text files [1]",
    )
    parser.add_argument("-o", "--outdir", type=str, default="", help="Directory to write to")
    parser.add_argument(
        "-z",
        "--procs",
        type=int,
        default=1,
        help="Number of parallel processes to run this job with",
    )
    parser.add_argument(
        "--vrna",
        type=str,
        default="",
        help="Append path to vrna RNA module to sys.path",
    )
    parser.add_argument(
        "-g",
        "--genes",
        type=str,
        default="",
        help="Genomic coordinates bed for genes, either standard bed format or AnnotateBed.pl format",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version",
    )
    parser.add_argument(
        "--loglevel",
        type=str,
        default="WARNING",
        choices=["WARNING", "ERROR", "INFO", "DEBUG"],
        help="Set log level",
    )
    parser.add_argument("--logdir", type=str, default="LOGS", help="Set log directory")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def parseargs_collectpl():
    parser = argparse.ArgumentParser(
        description="Calculate the regions with highest accessibility diff for given Sequence Pattern"
    )
    parser.add_argument(
        "-p",
        "--pattern",
        type=str,
        default="240,60",
        help="Pattern for files and window, e.g. Seq1_30,250",
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=1.0,
        help="Cutoff for the definition of pairedness, if set to > 0 it will select only constraint regions with mean raw (unconstraint) probability of being unpaired <= cutoff for further processing(default: 1.0)",
    )
    parser.add_argument(
        "-b",
        "--border",
        type=float,
        default=0.0,
        help="Cutoff for the minimum change between unconstraint and constraint structure, regions below this cutoff will not be further evaluated.",
    )
    parser.add_argument(
        "-u",
        "--ulimit",
        type=int,
        default=1,
        help="Stretch of nucleotides used during plfold run (-u option)",
    )
    parser.add_argument(
        "-r",
        "--unconstraint",
        type=str,
        default="raw",
        help="Name for unconstraint provided at ConstraintPLFold -r",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        type=int,
        default=37,
        help="Temperature for structure prediction",
    )
    parser.add_argument("-o", "--outdir", type=str, default="", help="Directory to write to")
    parser.add_argument("-d", "--dir", type=str, default="", help="Directory to read from")
    parser.add_argument(
        "-g",
        "--genes",
        type=str,
        help="Genomic coordinates bed for genes, either standard bed format or AnnotateBed.pl format",
    )
    parser.add_argument(
        "-z",
        "--procs",
        type=int,
        default=1,
        help="Number of parallel processes to run this job with",
    )
    parser.add_argument(
        "--loglevel",
        type=str,
        default="WARNING",
        choices=["WARNING", "ERROR", "INFO", "DEBUG"],
        help="Set log level",
    )
    parser.add_argument("--logdir", type=str, default="LOGS", help="Set log directory")
    parser.add_argument(
        "-w",
        "--padding",
        type=int,
        default=1,
        help="Padding around constraint that will be excluded from report, default is 1, "
        "so directly overlapping effects will be ignored",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def parseargs_browser():
    parser = argparse.ArgumentParser(description="Generates BigWig files for browsing from raw or constraint files")
    parser.add_argument(
        "-p",
        "--pattern",
        type=str,
        default="240,60",
        help="Pattern for and window and span, e.g. 30,250. Window can contain other strings for filtering, e.g. Seq1_30",
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=1.0,
        help="Cutoff for the definition of pairedness, if set to < 1 it will select only constraint regions with mean raw (unconstraint) probability of being unpaired <= cutoff for further processing(default: 1.0)",
    )
    parser.add_argument(
        "-b",
        "--border",
        type=float,
        default=0.0,
        help="Cutoff for the minimum change between unconstraint and constraint structure, regions below this cutoff will not be further evaluated.",
    )
    parser.add_argument(
        "-u",
        "--ulimit",
        type=int,
        default=1,
        help="Stretch of nucleotides used during plfold run (-u option)",
    )
    parser.add_argument(
        "-r",
        "--unconstraint",
        type=str,
        default=None,
        help="Name for unconstraint provided at ConstraintPLFold -r",
    )
    parser.add_argument(
        "-n",
        "--unpaired",
        action="store_true",
        help="If unpaired files should be converted as well",
    )
    parser.add_argument(
        "-a",
        "--paired",
        action="store_true",
        help="If paired files should be converted as well",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        type=int,
        default=37,
        help="Temperature for structure prediction",
    )
    parser.add_argument("-o", "--outdir", type=str, default="", help="Directory to write to")
    parser.add_argument("-d", "--dir", type=str, default="", help="Directory to read from")
    parser.add_argument(
        "-g",
        "--genes",
        type=str,
        help="Genomic coordinates bed for genes in standard BED format",
    )
    parser.add_argument(
        "-s",
        "--chromsizes",
        type=str,
        help="Chromosome sizes file",
    )
    parser.add_argument(
        "-z",
        "--procs",
        type=int,
        default=1,
        help="Number of parallel processes to run this job with",
    )
    parser.add_argument(
        "--loglevel",
        type=str,
        default="WARNING",
        choices=["WARNING", "ERROR", "INFO", "DEBUG"],
        help="Set log level",
    )
    parser.add_argument("--logdir", type=str, default="LOGS", help="Set log directory")
    parser.add_argument(
        "-w",
        "--padding",
        type=int,
        default=1,
        help="Padding around constraint that will be excluded from report, default is 1, "
        "so directly overlapping effects will be ignored",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def parseargs_foldcons():
    parser = argparse.ArgumentParser(
        description="Calculate base pairing probs of given seqs or random seqs for given window size, span and region."
    )
    parser.add_argument("-s", "--sequence", type=str, help="Sequence to fold")
    parser.add_argument("-w", "--window", type=int, default=None, help="Size of window")
    parser.add_argument("-l", "--span", type=int, default=None, help="Maximum bp span")
    parser.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=-0.01,
        help="Only print prob greater cutoff",
    )
    parser.add_argument(
        "-x",
        "--constrain",
        type=str,
        default="sliding",
        help="Region to constrain, either sliding window (default) or region to constrain (e.g. 1-10) "
        "or path to file containing regions following the naming pattern $fastaID_constraints, "
        "if paired, the first entry of the file will become a fixed constraint and paired "
        "with all the others, e.g. Sequence1_constraints, "
        "choices = [off, sliding, file, paired,"
        " or simply 1-10,2-11 or 1-10;15-20,2-11:16-21 for paired]",
    )
    parser.add_argument(
        "-y",
        "--conslength",
        type=int,
        default=0,
        help="Length of region to constrain for slidingwindow",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        type=int,
        default=37,
        help="Temperature for structure prediction",
    )
    parser.add_argument(
        "--save",
        type=str,
        default="STDOUT",
        help="Save the output as gz file with that name",
    )
    parser.add_argument("-o", "--outdir", type=str, default="", help="Directory to write to")
    parser.add_argument(
        "-z",
        "--procs",
        type=int,
        default=1,
        help="Number of parallel processes to run this job with",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="",
        help="Helper var, only used if called from other prog where a pattern for files is defined",
    )
    parser.add_argument(
        "-g",
        "--genes",
        type=str,
        default="",
        help="Genomic coordinates bed for genes in standard BED format",
    )
    parser.add_argument(
        "--loglevel",
        type=str,
        default="WARNING",
        choices=["WARNING", "ERROR", "INFO", "DEBUG"],
        help="Set log level",
    )
    parser.add_argument("--logdir", type=str, default="LOGS", help="Set log directory")
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def parseargs_collect_window():
    parser = argparse.ArgumentParser(
        description="Collect Fold output and calculate ddG for given sequence pattern and constraints"
    )
    parser.add_argument(
        "-p",
        "--pattern",
        type=str,
        default="250,150",
        help="Pattern for files and window, e.g. Seq1_30,250",
    )
    parser.add_argument(
        "-b",
        "--border",
        type=str,
        default="-inf,inf",
        help="Cutoff for the minimum change between unconstraint and constraint structure, regions below this cutoff will not be returned as list of regions with most impact on structure.",
    )
    parser.add_argument("-o", "--outdir", type=str, default="", help="Directory to write to")
    parser.add_argument(
        "-g",
        "--genes",
        type=str,
        help="Genomic coordinates bed for genes, standard BED format",
    )
    parser.add_argument(
        "-z",
        "--procs",
        type=int,
        default=1,
        help="Number of parallel processed to run this job with, only important of no border is given and we need to fold",
    )
    parser.add_argument(
        "--loglevel",
        type=str,
        default="WARNING",
        choices=["WARNING", "ERROR", "INFO", "DEBUG"],
        help="Set log level",
    )
    parser.add_argument("--logdir", type=str, default="LOGS", help="Set log directory")
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def parseargs_collect_windowdiff():
    parser = argparse.ArgumentParser(
        description="Calculate the regions with highest accessibility diff for given Sequence Pattern"
    )
    parser.add_argument(
        "-w",
        "--window",
        type=str,
        default="250",
        help="Pattern for window size, comma separated for list, e.g. 120,140,180,200",
    )
    parser.add_argument(
        "-b",
        "--span",
        type=str,
        default="150",
        help="Pattern for basepair span, comma separated for list, e.g. 120,140,180,200",
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        type=str,
        default="-inf,inf",
        help="Cutoff for the minimum change between unconstraint and constraint structure,"
        " regions below this cutoff will not be returned as list of regions with most "
        "impact on structure.",
    )
    parser.add_argument(
        "-u",
        "--ulimit",
        type=int,
        default=1,
        help="Stretch of nucleotides used during plfold run (-u option)",
    )
    parser.add_argument("-o", "--outdir", type=str, default="", help="Directory to write to")
    parser.add_argument(
        "-g",
        "--genes",
        type=str,
        help="Genomic coordinates bed for genes, either standard bed format or AnnotateBed.pl format",
    )
    parser.add_argument(
        "-z",
        "--procs",
        type=int,
        default=1,
        help="Number of parallel processed to run this job with,"
        " only important of no border is given and we need to fold",
    )
    parser.add_argument(
        "--loglevel",
        type=str,
        default="WARNING",
        choices=["WARNING", "ERROR", "INFO", "DEBUG"],
        help="Set log level",
    )
    parser.add_argument("--logdir", type=str, default="LOGS", help="Set log directory")
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def visualiziation_parser():
    parser = argparse.ArgumentParser(description="Visualize RIssmed Output using bed files")
    parser.add_argument("-f", "--file", type=str, help="Path to the RIssmed BED file")
    parser.add_argument(
        "-t",
        "--tmp",
        action="store_true",
        help="Will produce the database and store it in a tmp file",
    )
    parser.add_argument(
        "-m",
        "--memory",
        action="store_true",
        help="Will store database in memory overrides tmp flag",
    )
    parser.add_argument(
        "-d",
        "--database",
        type=str,
        help="path to store the database file. If it already exists " "this file is used and -f is ignored",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


#
# RIssmedArgparsers.py ends here
