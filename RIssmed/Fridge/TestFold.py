#!\usr\env\bin python

# TestFold.py ---

# IMPORTS
import argparse
import pprint
from io import StringIO
import time
import math
import os
import gzip
import importlib
import multiprocessing
from multiprocessing import Manager
import traceback as tb

# load own modules
import sys


from Collection import *
from Randseq import createrandseq

# Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq

# numpy and matplolib and pyplot
import numpy as np
import matplotlib
from matplotlib import animation  # , rc
import matplotTweaks.pyplot as plt
from random import choices, choice, shuffle  # need this if tempprobing was choosen

# find ffmpeg executable
# import shutil
# plt.rcParams['animation.ffmpeg_path'] = shutil.which("ffmpeg")
# plt.rc('verbose', level='debug-annoying', fileo=sys.stdout)
# matplotTweaks.verbose.set_level("helpful")
# plt.rc('animation', html='html5')


def parseargs():
    parser = argparse.ArgumentParser(
        description="Calculate base pairing probs of given seqs or random seqs for given window size, span and region."
    )
    parser.add_argument("-s", "--sequence", type=str, help="Sequence to fold")
    parser.add_argument("-w", "--window", type=int, default=240, help="Size of window")
    parser.add_argument("-l", "--span", type=int, default=60, help="Length of bp span")
    parser.add_argument("-u", "--region", type=int, default=5, help="Length of region")
    parser.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=-0.01,
        help="Only print prob greater cutoff",
    )
    parser.add_argument(
        "-r",
        "--unconstrained",
        type=str,
        default="STDOUT",
        help="Print output of unconstrained folding to file with this name",
    )
    parser.add_argument(
        "-n",
        "--unpaired",
        type=str,
        default="STDOUT",
        help="Print output of unpaired folding to file with this name",
    )
    parser.add_argument(
        "-p",
        "--paired",
        type=str,
        default="STDOUT",
        help="Print output of paired folding to file with this name",
    )
    parser.add_argument("-e", "--length", type=int, default=100, help="Length of randseq")
    parser.add_argument(
        "-g",
        "--gc",
        type=int,
        default=0,
        help="GC content, needs to be %2==0 or will be rounded",
    )
    parser.add_argument("-b", "--number", type=int, default=1, help="Number of random seqs to generate")
    parser.add_argument(
        "-x",
        "--constrain",
        type=str,
        default="sliding",
        help="Region to constrain, either sliding window (default) or region to constrain (e.g. 1-10) or path to file containing regions following the naming pattern $fastaID_constraints, if paired, the first entry of the file will become a fixed constraint and paired with all the others, e.g. Sequence1_constraints, choices = [off,sliding,temperature, tempprobe, file, paired, or simply 1-10,2-11 or 1-10;15-20,2-11;16-21 for paired]",
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
        "--temprange",
        type=str,
        default="",
        help="Temperature range for structure prediction (e.g. 37-60)",
    )
    parser.add_argument("-a", "--alphabet", type=str, default="AUCG", help="alphabet for random seqs")
    parser.add_argument(
        "--plot",
        type=str,
        default="0",
        choices=["0", "svg", "png"],
        help="Create image of the (un-)constraint sequence, you can select the file format here (svg,png). These images can later on be animated with ImageMagick like `convert -delay 120 -loop 0 *.svg animated.gif`.",
    )
    parser.add_argument("--save", type=int, default=1, help="Save the output as gz files")
    parser.add_argument("-o", "--outdir", type=str, default="", help="Directory to write to")
    parser.add_argument(
        "-z",
        "--procs",
        type=int,
        default=1,
        help="Number of parallel processed to run this job with",
    )
    parser.add_argument(
        "--vrna",
        type=str,
        default="",
        help="Append path to vrna RNA module to sys.path",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="",
        help="Helper var, only used if called from other prog where a pattern for files is defined",
    )
    parser.add_argument(
        "-v",
        "--verbosity",
        type=int,
        default=0,
        choices=[0, 1],
        help="increase output verbosity",
    )

    return parser.parse_args()


def fold(
    sequence,
    window,
    span,
    region,
    unconstrained,
    unpaired,
    paired,
    length,
    gc,
    number,
    constrain,
    conslength,
    alphabet,
    plot,
    save,
    procs,
    vrna,
    temprange,
    outdir,
    verbosity=False,
    pattern=None,
    cutoff=None,
):
    # set path for output
    if outdir:
        printlog(outdir)
        if not os.path.isabs(outdir):
            outdir = os.path.abspath(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    else:
        outdir = os.path.abspath(os.getcwd())
    # set path for VRNA lib
    if vrna:
        sys.path = [vrna] + sys.path
    else:
        sys.path = ["/scratch/fall/VRNA/243_alpha2/lib/python3.6/site-packages"] + sys.path
    try:
        global RNA
        RNA = importTweaks.import_module("RNA")
        globals().update(
            {n: getattr(RNA, n) for n in RNA.__all__}
            if hasattr(RNA, "__all__")
            else {k: v for (k, v) in RNA.__dict__.items() if not k.startswith("_")}
        )
    except ImportError:
        print("Error:", err)

    if plot == "0" and not save:
        raise ValueError(
            "Neither plot nor save are active, this script will take a long time and produce nothing, please activate at least one of the two!"
        )

    if sequence == "random":
        rand = "\n".join(createrandseq(length, gc, number, alphabet))
        seq = StringIO(rand)
        o = gzip.open("Random.fa.gz", "wb")
        o.write(bytes(rand, encoding="UTF-8"))
        o.close()

    elif os.path.isfile(sequence):
        if ".gz" in sequence:
            seq = gzip.open(sequence, "rt")
        else:
            seq = open(sequence, "rt")
    else:
        header = ">Seq1:default:nochrom:(.)"
        s = sequence
        seq = StringIO("{header}\n{s}".format(header=header, s=s))

    for fa in SeqIO.parse(seq, "fasta"):
        try:
            goi, chrom = fa.id.split(":")[::2]
            strand = str(fa.id.split(":")[3].split("(")[1][0])
        except:
            printlog("Fasta header is not in expected format, you will loose information on strand and chromosome")
        try:
            goi = fa.id
            chrom, strand = ["na", "na"]
        except:
            printlog("Could not assign any value from fasta header, please check your fasta files")

        if pattern and pattern not in goi:
            next
        else:
            printlog("Working on " + goi)
            ##prepare plots
            # 		if plot != '0':
            manager = Manager()
            animations = manager.list()
            xvals = []
            for y in range(1, len(fa.seq) + 1):
                xvals.append(y)
            xs = np.array(xvals)
            xvals = []
            # set kT for nrg2prob and vice versa calcs
            kT = 0.61632077549999997
            # define data structures
            # 			data = { 'bpp': [], 'up': [] }
            data = {"up": []}
            an = [None]
            # We check if we need to fold the whole seq or just a region around the constraints
            if constrain == "sliding" or constrain == "temperature":  # here we fold the whole seq
                # if not already available, run raw fold for whole sequence as we have a sliding window constraint
                check = (
                    str(goi) + "_" + str(chrom) + "_" + str(strand) + "_" + unconstrained + "_" + str(window) + ".gz"
                )
                if not (os.path.isfile(check)):
                    # set model details
                    md = RNA.md()
                    md.max_bp_span = span
                    md.window_size = window
                    # create new fold_compound object
                    fc = RNA.fold_compound(str(fa.seq), md, RNA.OPTION_WINDOW)
                    # call prop window calculation
                    fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data)
                    # 			fc.probs_window(region, RNA.PROBS_WINDOW_BPP, bpp_callback, data)
                    # print plfold output
                    if save:
                        write_unconstraint(fa, unconstrained, data, int(region), str(window), outdir)
                else:
                    # The hard work of folding this has already been done, so we just read in the results
                    printlog(
                        "Found "
                        + str(goi + "_" + unconstrained + "_" + str(window) + ".gz")
                        + ", will read in data and not recalculate"
                    )
                    data["up"] = read_precalc_fold(data["up"], check, fa)

                # convert to array for fast diff	calc
                an = up_to_array(data["up"], int(region), len(fa.seq))

            if constrain == "temperature":
                printlog("Calculating probs for temperature constraint " + temprange)
                # Create process pool with processes
                num_processes = procs or 1
                pool = multiprocessing.Pool(processes=num_processes)
                processes = []

                ts, te = map(int, temprange.split("-"))

                for temp in range(ts, te + 1):
                    # Create the process, and connect it to the worker function
                    try:
                        constrain_temp(
                            fa,
                            temp,
                            window,
                            span,
                            region,
                            an,
                            animations,
                            xs,
                            save,
                            outdir,
                            plot,
                        )
                    except Exception:
                        exc_type, exc_value, exc_tb = sys.exc_info()
                        tbe = tb.TracebackException(
                            exc_type,
                            exc_value,
                            exc_tb,
                        )
                        with open("error", "a") as h:
                            print("".join(tbe.format()), file=h)

            elif constrain == "tempprobe":
                # Create process pool with processes
                num_processes = procs or 1
                pool = multiprocessing.Pool(processes=num_processes)
                processes = []
                probeargs = temprange.split(",")
                printlog("Calculating cutoff probs for temperature constraint " + probeargs[0])
                ts, te = map(int, probeargs[0].split("-"))
                if (
                    len(fa.seq) > window * 4
                ):  # If the sequence is very large and we only need a proxy of changes to get a border for the cutoff of CalcConsDiffs we randomly select 10 regions of size window of the sequence for folding
                    selection = ""
                    if probeargs[1]:  # if we have constraints we choose windows around the constraints for tempprobing
                        for const in probeargs[1:]:
                            conss, conse = map(int, const.split("-"))
                            strt = conss - window * 2
                            endt = conse + window * 2 + 1
                            if strt < 0:
                                strt = 0
                            if endt > len(fa.seq):
                                endt = len(fa.seq)
                            selection = selection + str(fa.seq[strt:endt])
                    else:
                        for chosen in range(1, 11):  # we randomly choose 10 windows from the sequence
                            items = range(window, len(fa.seq) - window - 1)
                            sel = choice(items)
                            selection = selection + (fa.seq[sel : sel + window])
                    fa.seq = Seq(str(selection))
                # we check if we already have the raw seq folded
                check = (
                    str(goi) + "_" + str(chrom) + "_" + str(strand) + "_" + unconstrained + "_" + str(window) + ".gz"
                )
                if os.path.isfile(check):
                    # The hard work of folding this has already been done, so we just read in the results
                    printlog("Found " + check + ", will read in data and not recalculate")
                    data["up"] = read_precalc_fold(data["up"], check, fa)
                    # convert to array for fast diff	calc
                    an = up_to_array(data["up"], int(region), len(fa.seq))

                if an[0] == None or len(an) > len(
                    fa.seq
                ):  # This means that the prob info was loaded from file or the raw sequence has not been folded yet, but we need a subsequence for the cutoff calculation, so we fold the subsequence
                    printlog("Recalculating at default temp with subseq")
                    # set model details
                    md = RNA.md()
                    md.max_bp_span = span
                    md.window_size = window

                    # create new fold_compound object
                    fc = RNA.fold_compound(str(fa.seq), md, RNA.OPTION_WINDOW)
                    # call prop window calculation
                    fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data)
                    an = up_to_array(data["up"], int(region), len(fa.seq))

                for temp in range(ts, te + 1):
                    # 					print('Recalculating at temp '+str(temp))
                    # Create the process, and connect it to the worker function
                    try:
                        constrain_temp(
                            fa,
                            temp,
                            window,
                            span,
                            region,
                            an,
                            animations,
                            xs,
                            save,
                            outdir,
                            plot,
                        )
                    except Exception:
                        exc_type, exc_value, exc_tb = sys.exc_info()
                        tbe = tb.TracebackException(
                            exc_type,
                            exc_value,
                            exc_tb,
                        )
                        with open("error", "a") as h:
                            print("".join(tbe.format()), file=h)

            else:
                # Create process pool with processes
                num_processes = procs or 1
                pool = multiprocessing.Pool(processes=num_processes)
                conslist = []

                if os.path.isfile(constrain):
                    printlog("Calculating probs for constraint from file " + constrain)
                    with gzip.open(constrain, "rt") as o:
                        for line in o:
                            conslist.append(line.rstrip())
                elif constrain == "file" or constrain == "paired":
                    printlog("Calculating probs for constraint from file " + str(goi + "_constraints"))
                    with open(goi + "_constraints", "rt") as o:
                        for line in o:
                            conslist.append(line.rstrip())
                elif constrain == "none":
                    conslist = ["NOCONS"]
                elif constrain == "sliding":
                    for start in range(1, len(fa.seq) - conslength + 2):
                        end = start + conslength - 1
                        conslist.append(str(start) + "-" + str(end))
                else:
                    printlog("Calculating probs for constraint " + constrain)
                    conslist = constrain.split(",")

                for entry in conslist:
                    if entry == "NOCONS":  # in case we just want to fold the sequence without constraints at all
                        md = RNA.md()
                        md.max_bp_span = span
                        md.window_size = window

                        # create new fold_compound object
                        fc = RNA.fold_compound(str(fa.seq), md, RNA.OPTION_WINDOW)
                        # call prop window calculation
                        data = {"up": []}
                        fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data)
                        write_unconstraint(fa, unconstrained, data, int(region), str(window), outdir)

                    else:
                        # we now have a list of constraints and for the raw seq comparison we only need to fold windows around these constraints
                        # In case we want to constrain pairwise
                        fstart, fend = [None, None]
                        data = {"up": []}
                        if constrain == "paired" or ":" in entry:
                            [fstart, fend], [start, end] = [
                                [int(x) for x in cn.split("-", 1)] for cn in entry.split(":", 1)
                            ]
                            tostart = fstart - 4 * window
                            if tostart < 0:
                                tostart = 0
                            toend = fend + 4 * window + 1
                            if toend > len(fa.seq):
                                toend = len(fa.seq)
                            cons = str(fstart) + "-" + str(fend) + ":" + str(start) + "-" + str(end)
                        else:
                            start, end = map(int, entry.split("-", 1))
                            # for all constraints we now extract subsequences to compare against
                            # we no longer fold the whole raw sequence but only the constraint region +- window size
                            tostart = start - 4 * window
                            if tostart < 0:
                                tostart = 0
                            toend = end + 4 * window + 1
                            if toend > len(fa.seq):
                                toend = len(fa.seq)
                            cons = str(start) + "-" + str(end)

                        if checkexisting(fa, paired, unpaired, cons, region, window, outdir):
                            # 							printlog(cons + ' EXIST')
                            continue

                        printlog("Calculating constraint\t" + entry)
                        seqtofold = str(fa.seq[tostart:toend])  ###TEST
                        const = np.array([start, end])

                        # 						print('Constraint: ' + str(entry) + '\tSeqlength: ' + str(len(fa.seq)) + '\tFoldlength: ' + str(len(seqtofold)) + ' : ' + str(tostart) + " - " + str(toend))
                        # Now we fold this region to get the raw probs
                        # set model details
                        md = RNA.md()
                        md.max_bp_span = span
                        md.window_size = window

                        # create new fold_compound object
                        fc = RNA.fold_compound(str(seqtofold), md, RNA.OPTION_WINDOW)
                        # In case we have paired constraint we want to set the first constraint as default
                        # call prop window calculation
                        data_t = {"up": []}
                        fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_t)
                        # for easier handling we now extend the array before and after the constraint region to the total length of the raw sequence
                        for i in range(0, tostart):
                            data["up"].append([None, 0.0])
                        data["up"].extend(data_t["up"][:])
                        data_t = None  # remove unneeded list
                        for i in range(toend, len(fa.seq)):
                            data["up"].append([None, 0.0])

                        an = up_to_array(data["up"], int(region), len(fa.seq))

                        if fstart is not None and fend is not None:
                            printlog("Constraining to " + str(fstart) + " and " + str(fend))
                            try:
                                goi, chrom = fa.id.split(":")[::2]
                                strand = str(fa.id.split(":")[3].split("(")[1][0])
                            except:
                                goi = fa.id
                                chrom, strand = ["na", "na"]
                            if not os.path.exists(
                                "StruCons_"
                                + goi
                                + "_"
                                + chrom
                                + "_"
                                + strand
                                + "_"
                                + cons
                                + "_"
                                + paired
                                + "_"
                                + str(window)
                                + ".gz"
                            ) and not os.path.exists(
                                "StruCons_"
                                + goi
                                + "_"
                                + chrom
                                + "_"
                                + strand
                                + "_"
                                + cons
                                + "_"
                                + unpaired
                                + "_"
                                + str(window)
                                + ".gz"
                            ):
                                try:
                                    constrain_seq(
                                        fa,
                                        fstart,
                                        fend,
                                        conslength,
                                        const,
                                        str(fstart) + "-" + str(fend),
                                        window,
                                        span,
                                        region,
                                        an,
                                        animations,
                                        xs,
                                        paired,
                                        unpaired,
                                        save,
                                        outdir,
                                        plot,
                                        data,
                                    )
                                except Exception:
                                    exc_type, exc_value, exc_tb = sys.exc_info()
                                    tbe = tb.TracebackException(
                                        exc_type,
                                        exc_value,
                                        exc_tb,
                                    )
                                    with open("error", "a") as h:
                                        print("".join(tbe.format()), file=h)
                            try:
                                constrain_seq_paired(
                                    fa,
                                    fstart,
                                    fend,
                                    start,
                                    end,
                                    conslength,
                                    const,
                                    cons,
                                    window,
                                    span,
                                    region,
                                    an,
                                    animations,
                                    xs,
                                    paired,
                                    unpaired,
                                    save,
                                    outdir,
                                    plot,
                                    data,
                                )
                            except Exception:
                                exc_type, exc_value, exc_tb = sys.exc_info()
                                tbe = tb.TracebackException(
                                    exc_type,
                                    exc_value,
                                    exc_tb,
                                )
                                with open("error", "a") as h:
                                    print("".join(tbe.format()), file=h)
                        else:
                            try:
                                constrain_seq(
                                    fa,
                                    start,
                                    end,
                                    conslength,
                                    const,
                                    cons,
                                    window,
                                    span,
                                    region,
                                    an,
                                    animations,
                                    xs,
                                    paired,
                                    unpaired,
                                    save,
                                    outdir,
                                    plot,
                                    data,
                                )
                            except Exception:
                                exc_type, exc_value, exc_tb = sys.exc_info()
                                tbe = tb.TracebackException(
                                    exc_type,
                                    exc_value,
                                    exc_tb,
                                )
                                with open("error", "a") as h:
                                    print("".join(tbe.format()), file=h)

    printlog("DONE: output in: " + str(outdir))


##### Functions #####
def constrain_seq(
    fa,
    start,
    end,
    conslength,
    const,
    cons,
    window,
    span,
    region,
    an,
    animations,
    xs,
    paired,
    unpaired,
    save,
    outdir,
    plot,
    data,
):
    #   DEBUGGING
    # 	pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)

    try:
        goi, chrom = fa.id.split(":")[::2]
        strand = str(fa.id.split(":")[3].split("(")[1][0])
    except:
        goi = fa.id
        chrom, strand = ["na", "na"]

    if os.path.exists(
        "StruCons_" + goi + "_" + chrom + "_" + strand + "_" + cons + "_" + paired + "_" + str(window) + ".gz"
    ) and os.path.exists(
        "StruCons_" + goi + "_" + chrom + "_" + strand + "_" + cons + "_" + unpaired + "_" + str(window) + ".gz"
    ):
        return

    # we no longer fold the whole sequence but only the constraint region +- window size
    tostart = start - 4 * window
    if tostart < 0:
        tostart = 0
    toend = end + 4 * window + 1
    if toend > len(fa.seq):
        toend = len(fa.seq)
    seqtofold = str(fa.seq[tostart:toend])  ###TEST

    # refresh model details
    md = RNA.md()
    md.max_bp_span = span
    md.window_size = window

    # create new fold_compound objects
    fc_p = RNA.fold_compound(seqtofold, md, RNA.OPTION_WINDOW)
    fc_u = RNA.fold_compound(seqtofold, md, RNA.OPTION_WINDOW)

    # enforce paired
    for x in range(start - tostart, end - tostart + 1):
        fc_p.hc_add_bp_nonspecific(
            x, 0
        )  # 0 means without direction  ( $ d < 0 $: pairs upstream, $ d > 0 $: pairs downstream, $ d == 0 $: no direction)

    # enforce unpaired
    for x in range(start - tostart, end - tostart + 1):
        fc_u.hc_add_up(x)

    # new data struct
    data_pn = {"up": []}
    data_un = {"up": []}

    fc_p.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_pn)
    fc_u.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_un)

    # 	with open('bla','a') as h:
    # 		print(str(tostart)+'\t'+str(toend)+'\t'+str(len(data['up']))+cons+'\n'+str(data['up'][tostart:toend])+'\n'+str(data_pn['up'][:window+window+end-start+1])+'\n'+str(data_un['up'][:window+window+end-start+1]), file=h)

    # we now fill the list with values from the unconstrained sequence to get the same length
    data_p = {"up": []}
    data_u = {"up": []}

    data_p["up"].extend(data["up"][0:tostart])
    data_p["up"].extend(data_pn["up"][:])
    data_pn = None  # remove unneeded list
    data_p["up"].extend(data["up"][toend:])

    data_u["up"].extend(data["up"][0:tostart])
    data_u["up"].extend(data_un["up"][:])
    data_un = None  # remove unneeded list
    data_u["up"].extend(data["up"][toend:])

    # 	with open('bla2','a') as h:
    # 		print(str(data['up'][tostart-1:tostart+9])+cons+'\n'+str(data_p['up'][tostart-1:tostart+9])+'\n'+str(data_u['up'][tostart-1:tostart+9])+'\n'+str(len(data['up']))+'\t'+str(len(data_p['up']))+'\t'+str(len(data_u['up'])), file=h)

    au = up_to_array(data_u["up"], region, len(fa.seq))
    ap = up_to_array(data_p["up"], region, len(fa.seq))

    if not np.array_equal(an, au):
        diff_nu = an - au
    else:
        printlog("No influence on Structure with unpaired constraint at " + cons)
        diff_nu = None
    if not np.array_equal(an, ap):
        diff_np = an - ap
    else:
        printlog("No influence on Structure with paired constraint at " + cons)
        diff_np = None

    if plot == "svg" or plot == "png":
        # 		printlog('PLOTTING '+str(plot))
        plot_data(fa, an, au, ap, const, xs, cons, plot, outdir)

    if save:
        # 		printlog('SAVING')
        write_constraint(
            fa,
            paired,
            unpaired,
            data_u,
            data_p,
            cons,
            region,
            diff_nu,
            diff_np,
            str(window),
            outdir,
        )


def constrain_seq_paired(
    fa,
    fstart,
    fend,
    start,
    end,
    conslength,
    const,
    cons,
    window,
    span,
    region,
    an,
    animations,
    xs,
    paired,
    unpaired,
    save,
    outdir,
    plot,
    data,
):
    #   DEBUGGING
    # 	pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)

    # we no longer fold the whole sequence but only the constraint region +- window size
    tostart = start - 4 * window
    if tostart < 0:
        tostart = 0
    toend = end + 4 * window + 1
    if toend > len(fa.seq):
        toend = len(fa.seq)
    seqtofold = str(fa.seq[tostart:toend])  ###TEST

    # refresh model details
    md = RNA.md()
    md.max_bp_span = span
    md.window_size = window

    # create new fold_compound objects
    fc_p = RNA.fold_compound(seqtofold, md, RNA.OPTION_WINDOW)
    fc_u = RNA.fold_compound(seqtofold, md, RNA.OPTION_WINDOW)

    # enforce paired
    for x in range(start - tostart, end - tostart + 1):
        fc_p.hc_add_bp_nonspecific(
            x, 0
        )  # 0 means without direction  ( $ d < 0 $: pairs upstream, $ d > 0 $: pairs downstream, $ d == 0 $: no direction)
        # enforce paired
    for x in range(fstart - tostart, fend - tostart + 1):
        fc_p.hc_add_bp_nonspecific(x, 0)  # 0 means without direction  ( $ d < 0 $: pairs upstream, $ d >

    # enforce unpaired
    for x in range(start - tostart, end - tostart + 1):
        fc_u.hc_add_up(x)
    # enforce unpaired
    for x in range(fstart - tostart, fend - tostart + 1):
        fc_u.hc_add_up(x)

    # new data struct
    data_pn = {"up": []}
    data_un = {"up": []}

    fc_p.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_pn)
    fc_u.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_un)

    # 	with open('bla','a') as h:
    # 		print(str(tostart)+'\t'+str(toend)+'\t'+str(len(data['up']))+'\n'+str(data['up'][tostart:toend])+'\n'+str(data_pn['up'][:window+window+end-start+1])+'\n'+str(data_un['up'][:window+window+end-start+1]), file=h)

    # we now fill the list with values from the unconstrained sequence to get the same length
    data_p = {"up": []}
    data_u = {"up": []}

    data_p["up"].extend(data["up"][0:tostart])
    data_p["up"].extend(data_pn["up"][:])
    data_pn = None  # remove unneeded list
    data_p["up"].extend(data["up"][toend:])

    data_u["up"].extend(data["up"][0:tostart])
    data_u["up"].extend(data_un["up"][:])
    data_un = None  # remove unneeded list
    data_u["up"].extend(data["up"][toend:])

    # 	with open('bla2','a') as h:
    # 		print(str(data['up'][tostart-1:tostart+20])+'\n'+str(data_p['up'][tostart-1:tostart+20])+'\n'+str(data_u['up'][tostart-1:tostart+20])+'\n'+str(len(data['up']))+'\t'+str(len(data_p['up']))+'\t'+str(len(data_u['up'])), file=h)

    au = up_to_array(data_u["up"], region, len(fa.seq))
    ap = up_to_array(data_p["up"], region, len(fa.seq))

    if not np.array_equal(an, au):
        diff_nu = an - au
    else:
        printlog("No influence on Structure with unpaired constraint at " + cons)
        diff_nu = None
    if not np.array_equal(an, ap):
        diff_np = an - ap
    else:
        printlog("No influence on Structure with paired constraint at " + cons)
        diff_np = None

    if plot == "svg" or plot == "png":
        # 		printlog('PLOTTING '+str(plot))
        plot_data(fa, an, au, ap, const, xs, cons, plot, outdir)

    if save:
        # 		printlog('SAVING')
        write_constraint(
            fa,
            paired,
            unpaired,
            data_u,
            data_p,
            cons,
            region,
            diff_nu,
            diff_np,
            str(window),
            outdir,
        )


def constrain_temp(fa, temp, window, span, region, an, animations, xs, save, outdir, plot):
    # 	print('FOLDING ' + str(fa.seq) + ' at temp ' + str(temp))
    # refresh model details
    md = RNA.md()
    md.max_bp_span = span
    md.window_size = window
    # set temperature
    md.temperature = temp
    # create new fold_compound objects
    fc_t = RNA.fold_compound(str(fa.seq), md, RNA.OPTION_WINDOW)
    # new data struct
    data_t = {"up": []}

    fc_t.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data_t)

    at = up_to_array(data_t["up"], region, len(fa.seq))

    diff_nt = an - at

    #####Set temp, and rewrite save and plot for temp
    if save:
        # 		print('SAVING')
        write_temp(fa, str(temp), data_t, region, diff_nt, str(window), outdir)

    if plot == "svg" or plot == "png":
        # 		print('PLOTTING '+str(plot))
        plot_temp(fa, at, temp, xs, plot, outdir)


def plot_data(fa, raw, consu, consp, const, xs, cons, saveas, outdir):

    anime = []
    # define xs for constraint line
    consl = []
    for x in const:
        consl.append(1.25)
    width = 16 / 100 * len(fa.seq)
    height = 9
    fig = plt.figure(figsize=(width, height), dpi=80)
    ax1 = fig.add_subplot(111)

    ax2 = ax1.twiny()
    # 	line, = ax.plot([], [], lw=2)
    plt.title(
        "Blue-- = unconstrained, Green-. = Unpaired, Red = Paired, Gray = Constraint",
        y=1.075,
    )
    ax1.set_ylabel("Prob unpaired")
    ax1.set_xlabel("Nucleotides")
    # 	plt.xticks(range(0,len(fa.seq)+1),(' '+fa.seq),size='small')
    # add lines to plot
    ax1.plot(xs, raw, "b-", xs, consu, "g-", xs, consp, "r-", const, consl, "k-")
    ax1.set_xlim(0, len(fa.seq) + 1)
    ax2.set_xlim(ax1.get_xlim())
    ax1.set_xticks(range(0, len(fa.seq) + 1))
    ax1.set_xticklabels((" " + fa.seq), ha="right")
    # 	ax2.set_xlabel(r"Modified x-axis: $1/(1+X)$")
    ax2.set_xticks(range(1, len(fa.seq) + 1))
    ax2.set_xticklabels(range(1, len(fa.seq) + 1), rotation=45, ha="right")
    # We change the fontsize of minor ticks label
    ax1.tick_params(axis="both", which="major", labelsize=8)
    ax1.tick_params(axis="both", which="minor", labelsize=4)
    ax2.tick_params(axis="both", which="major", labelsize=5)
    ax2.tick_params(axis="both", which="minor", labelsize=3)
    goi, chrom = fa.id.split(":")[::2]
    strand = str(fa.id.split(":")[3].split("(")[1][0])
    fig.savefig("StruCons_" + goi + "_" + cons + "." + saveas)
    plt.close()


# 	anime.append(plt.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-'))
# 	return anime


def plot_temp(fa, raw, temp, xs, saveas, outdir):

    anime = []
    # define xs for constraint line
    width = 16 / 100 * len(fa.seq)
    height = 9
    fig = plt.figure(figsize=(width, height), dpi=80)
    ax1 = fig.add_subplot(111)
    plt.title("Blue-- = " + temp + " degree", y=1.075)
    ax1.set_ylabel("Prob unpaired")
    ax1.set_xlabel("Nucleotides")
    # add lines to plot
    ax1.plot(xs, raw, "b-")
    ax1.set_xlim(0, len(fa.seq) + 1)
    ax1.set_xticks(range(0, len(fa.seq) + 1))
    ax1.set_xticklabels((" " + fa.seq))
    # We change the fontsize of minor ticks label
    ax1.tick_params(axis="both", which="major", labelsize=8)
    ax1.tick_params(axis="both", which="minor", labelsize=4)
    goi, chrom = fa.id.split(":")[::2]
    strand = str(fa.id.split(":")[3].split("(")[1][0])
    fig.savefig("TempCons_" + goi + "_" + temp + "." + saveas)
    plt.close()


# 	anime.append(plt.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-'))
# 	return anime


def bpp_callback(v, v_size, i, maxsize, what, data):
    if what & RNA.PROBS_WINDOW_BPP:
        data["bpp"].extend([{"i": i, "j": j, "p": p} for j, p in enumerate(v) if (p is not None)])  # and (p >= 0.01)])


def up_callback(v, v_size, i, maxsize, what, data):
    if what & RNA.PROBS_WINDOW_UP:
        #        data['up'].extend([{ 'i': i, 'up': v}])
        data["up"].extend([v])


def print_up(data=None, seqlength=None, region=None):
    # 	pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    # 	pp.pprint(data)
    try:
        ups = ""
        x = int(region)
        print("Printing region " + str(region) + " of data " + str(data[0]))
        for i in range(int(seqlength)):
            if data[i][x] is None or data[i][x] is np.nan or data[i][x] is "nan":
                data[i][x] = "NA"
            else:
                data[i][x] = round(data[i][x], 7)
            ups += str(i + 1) + "\t" + str(data[i][x]) + "\n"
        return ups
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def print_region_up(data=None, seqlength=None, region=None):
    # 	pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    # 	pp.pprint(data)
    try:
        ups = ""
        for i in range(int(seqlength)):
            for x in range(1, region + 1):
                if data[i][x] is None or data[i][x] is np.nan:
                    data[i][x] = "NA"
                else:
                    data[i][x] = round(data[i][x], 7)

            ups += str(i + 1) + "\t" + "\t".join(map(str, data[i][1 : region + 1])) + "\n"
        return ups
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def up_to_array(data=None, region=None, seqlength=None):
    # 	pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    # 	pp.pprint(data[165553:165588])
    try:
        entries = []
        if not seqlength:
            seqlength = len(data)

        for i in range(seqlength):
            if data[i][region] is None or data[i][region] is "NA" or data[i][region] is "nan":
                data[i][region] = np.nan
                entries.append(data[i][region])
            else:
                entries.append(round(data[i][region], 7))

        return np.array(entries)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def npprint(a, o=None):  # , format_string ='{0:.2f}'):
    try:
        out = ""
        it = np.nditer(a, flags=["f_index"])
        while not it.finished:
            out += "%d\t%0.7f" % (it.index + 1, it[0]) + "\n"
            it.iternext()
        if o:
            o.write(bytes(out, encoding="UTF-8"))
        else:
            print(out)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def write_unconstraint(fa, unconstrained, data, region, window, outdir):
    # 	pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    # 	pp.pprint(data['up'][1:10])
    try:
        goi, chrom = fa.id.split(":")[::2]
        strand = str(fa.id.split(":")[3].split("(")[1][0])
    except:
        goi = fa.id
        chrom, strand = ["na", "na"]

    try:
        if unconstrained != "STDOUT":
            bakdir = os.path.abspath(os.getcwd())
            os.chdir(outdir)
            if not os.path.exists(str(goi + "_" + chrom + "_" + strand + "_" + unconstrained + "_" + window + ".gz")):
                o = gzip.open(
                    goi + "_" + chrom + "_" + strand + "_" + unconstrained + "_" + window + ".gz",
                    "wb",
                )
                o.write(bytes(print_up(data["up"], len(fa.seq), region), encoding="UTF-8"))
            os.chdir(bakdir)
        else:
            print(print_up(data, len(fa.seq), region))
        return 1
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def write_constraint(
    fa,
    paired,
    unpaired,
    data_u,
    data_p,
    constrain,
    region,
    diff_nu,
    diff_np,
    window,
    outdir,
):
    # 	pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    # 	pp.pprint(data_u['up'][1:10])

    try:
        goi, chrom = fa.id.split(":")[::2]
        strand = str(fa.id.split(":")[3].split("(")[1][0])
    except:
        goi = fa.id
        chrom, strand = ["na", "na"]

    # print outputs to file or STDERR
    try:
        if paired != "STDOUT":
            bakdir = os.path.abspath(os.getcwd())
            os.chdir(outdir)
            if not os.path.exists(
                "StruCons_" + goi + "_" + chrom + "_" + strand + "_" + constrain + "_" + paired + "_" + window + ".gz"
            ):
                o = gzip.open(
                    "StruCons_"
                    + goi
                    + "_"
                    + chrom
                    + "_"
                    + strand
                    + "_"
                    + constrain
                    + "_"
                    + paired
                    + "_"
                    + window
                    + ".gz",
                    "wb",
                )
                o.write(bytes(print_up(data_p["up"], len(fa.seq), region), encoding="UTF-8"))
            os.chdir(bakdir)
        else:
            print(print_up(data_p["up"], len(fa.seq), region))

        if unpaired != "STDOUT":
            bakdir = os.path.abspath(os.getcwd())
            os.chdir(outdir)
            if not os.path.exists(
                "StruCons_"
                + goi
                + "_"
                + chrom
                + "_"
                + strand
                + "_"
                + constrain
                + "_"
                + unpaired
                + "_"
                + window
                + ".gz"
            ):
                o = gzip.open(
                    "StruCons_"
                    + goi
                    + "_"
                    + chrom
                    + "_"
                    + strand
                    + "_"
                    + constrain
                    + "_"
                    + unpaired
                    + "_"
                    + window
                    + ".gz",
                    "wb",
                )
                o.write(bytes(print_up(data_u["up"], len(fa.seq), region), encoding="UTF-8"))
            os.chdir(bakdir)
        else:
            print(print_up(data_u["up"], len(fa.seq), region))

        if diff_nu.any():
            if unpaired != "STDOUT":
                bakdir = os.path.abspath(os.getcwd())
                os.chdir(outdir)
                if not os.path.exists(
                    "StruCons_" + goi + "_" + chrom + "_" + strand + "_" + constrain + "_diffnu_" + window + ".gz"
                ):
                    o = gzip.open(
                        "StruCons_" + goi + "_" + chrom + "_" + strand + "_" + constrain + "_diffnu_" + window + ".gz",
                        "wb",
                    )
                    npprint(diff_nu, o)
                os.chdir(bakdir)
            else:
                npprint(diff_nu)

        if diff_np.any():
            if unpaired != "STDOUT":
                bakdir = os.path.abspath(os.getcwd())
                os.chdir(outdir)
                if not os.path.exists(
                    "StruCons_" + goi + "_" + chrom + "_" + strand + "_" + constrain + "_diffnp_" + window + ".gz"
                ):
                    o = gzip.open(
                        "StruCons_" + goi + "_" + chrom + "_" + strand + "_" + constrain + "_diffnp_" + window + ".gz",
                        "wb",
                    )
                    npprint(diff_np, o)
                os.chdir(bakdir)
            else:
                npprint(diff_np)
        return 1
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def write_temp(fa, temp, data, region, diff, window, outdir):
    # print outputs to file or STDERR
    try:
        goi, chrom = fa.id.split(":")[::2]
        strand = str(fa.id.split(":")[3].split("(")[1][0])
    except:
        goi = fa.id
        chrom, strand = ["na", "na"]

    try:
        bakdir = os.path.abspath(os.getcwd())
        os.chdir(outdir)

        if not os.path.exists("TempCons_" + goi + "_" + chrom + "_" + strand + "_" + temp + "_temp_" + window + ".gz"):
            o = gzip.open(
                "TempCons_" + goi + "_" + chrom + "_" + strand + "_" + temp + "_temp_" + window + ".gz",
                "wb",
            )
            o.write(bytes(print_up(data["up"], len(fa.seq), region), encoding="UTF-8"))
            o.close()
            if not os.path.exists(
                "TempCons_" + goi + "_" + chrom + "_" + strand + "_" + temp + "_difft_" + window + ".gz"
            ):
                o = gzip.open(
                    "TempCons_" + goi + "_" + chrom + "_" + strand + "_" + temp + "_difft_" + window + ".gz",
                    "wb",
                )
                npprint(diff, o)
                o.close()
        os.chdir(bakdir)
        return 1
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def read_precalc_fold(data, name, fa):
    try:
        for i in range(len(fa.seq)):
            data.append([])
            data[i] = []
        with gzip.open(name, "rt") as o:
            for line in o:
                cells = line.rstrip().split("\t")
                data[int(cells[0]) - 1].append([])
                data[int(cells[0]) - 1][0] = None
                for a in range(1, len(cells)):
                    data[int(cells[0]) - 1].append([])
                    data[int(cells[0]) - 1][a] = float(cells[a])
        return data
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def checkexisting(fa, paired, unpaired, cons, region, window, outdir):
    try:
        goi, chrom = fa.id.split(":")[::2]
        strand = str(fa.id.split(":")[3].split("(")[1][0])
    except:
        goi = fa.id
        chrom, strand = ["na", "na"]
    try:
        if os.path.exists(
            "StruCons_" + goi + "_" + chrom + "_" + strand + "_" + cons + "_" + paired + "_" + str(window) + ".gz"
        ) and os.path.exists(
            "StruCons_" + goi + "_" + chrom + "_" + strand + "_" + cons + "_" + unpaired + "_" + str(window) + ".gz"
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
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


if __name__ == "__main__":
    args = parseargs()
    fold(
        args.sequence,
        args.window,
        args.span,
        args.region,
        args.unconstrained,
        args.unpaired,
        args.paired,
        args.length,
        args.gc,
        args.number,
        args.constrain,
        args.conslength,
        args.alphabet,
        args.plot,
        args.save,
        args.procs,
        args.vrna,
        args.temprange,
        args.outdir,
        args.verbosity,
        args.pattern,
        args.cutoff,
    )
# ConstraintPLFold.py ends here

#
# TestPLFold.py ends here
