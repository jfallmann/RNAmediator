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
#find ffmpeg executable
#import shutil
#plt.rcParams['animation.ffmpeg_path'] = shutil.which("ffmpeg")
#plt.rc('verbose', level='debug-annoying', fileo=sys.stdout)
#matplotlib.verbose.set_level("helpful")
#plt.rc('animation', html='html5')
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
from RIssmed.RNAtweaks.logger import makelogdir, makelogfile, listener_process, listener_configurer, worker_configurer
# load own modules
from RIssmed.RNAtweaks.FileProcessor import *
from RIssmed.RNAtweaks.RNAtweaks import *
from lib.NPtweaks import *


log = logging.getLogger(__name__)  # use module name
scriptname = os.path.basename(__file__).replace('.py', '')


def fold(sequence, window, span, unconstraint, unpaired, paired, constrain, conslength, save, procs, outdir, run_settings, pattern=None, cutoff=None, queue=None, configurer=None, level=None):

    logid = scriptname+'.fold: '
    try:
        if queue and level:
            configurer(queue, level)

        # Create process pool with processes
        num_processes = procs or 1
        #with get_context("spawn").Pool(processes=num_processes-1, maxtasksperchild=1) as pool:
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
            log.debug(logid+str(conslist))

            if pattern and pattern not in goi:
                continue
            else:
                log.info(logid+'Working on ' + goi + "\t" + fa.id)

            # constraints
            for entry in conslist:
                log.debug(logid + 'ENTRY: '+str(entry))

                if not window:
                    window = len(fa.seq)
                if not span:
                    span = len(fa.seq)

                if len(fa.seq) < window:
                    log.warning('Sequence of '+goi+' to short, seqlenght '+str(len(fa.seq))+' with window size '+str(window))
                    continue

                if entry == 'NOCONS': # in case we just want to fold the sequence without constraints at all
                    gibbs_uc = [pool.apply_async(fold_unconstraint, args=(fa.seq), kwds={'queue':queue, 'configurer':configurer, 'level':level})]
                    return (gibbs_uc)

                else:
                    # we now have a list of constraints and for the raw seq comparison we only need to fold windows around these constraints
                    fstart, fend = [None,None]
                    if any(x in constrain for x in ['paired','Paired']) or ':' in entry:  # Not strand dependend, still genomic coords
                        if gstrand == '+' or gstrand == '.':
                            [fstart, fend], [start, end] = [[x - gs for x in get_location(cn)[:2]] for cn in entry.split(':',1)]
                        else:
                            [fstart, fend], [start, end] = [[ge - x for x in get_location(cn)[:2][::-1]] for cn in entry.split(':',1)]
                        cons = str(fstart)+'-'+str(fend)+':'+str(start)+'-'+str(end)
                        const = np.array([fstart, fend, start, end])
                        if start < 0 or fstart < 0 or end > len(fa.seq) or fend > len(fa.seq):
                            log.warning(logid+'Constraint out of sequence bounds! skipping! '+','.join(map(str,[goi,len(fa.seq),str(start)+'-'+str(end),str(fstart)+'-'+str(fend)])))
                            continue

                    else:
                        if gstrand == '+' or gstrand == '.':
                            start, end = [x - gs for x in get_location(entry)[:2]]
                        else:
                            start, end = [ge - x for x in get_location(entry)[:2][::-1]]

                        cons = str(start)+'-'+str(end)
                        log.debug(logid+str.join(' ',[goi,cons,gstrand]))
                        const = np.array([start, end])

                    pool.apply_async(constrain_seq, args=(fa, start, end, conslength, const, cons, window, span, unconstraint, paired, unpaired, save, outdir, genecoords), kwds={'queue':queue, 'configurer':configurer, 'level':level})

            pool.close()
            pool.join()

        log.info(logid+"DONE: output in: " + str(outdir))

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


##### Functions #####


def constrain_seq(fa, start, end, conslength, const, cons, window, span, unconstraint, paired, unpaired, save, outdir, genecoords, queue=None, configurer=None, level=None):
    #   DEBUGGING
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)

    logid = scriptname+'.constrain_seq: '
    try:
        if queue and level:
            configurer(queue, level)

        outlist = list()
        goi, chrom, strand = idfromfa(fa.id)

        if not window:
            window = len(fa.seq)
        if not span:
            span = len(fa.seq)

        dist, fstart, fend = [None, None, None]
        if ':' in cons:
            start, end, fstart, fend = const
            dist=abs(start-fend)
            if dist > window:
                return(log.warning(logid+'Window '+str(window)+' too small for constraint distance '+str(dist)))
        else:
            start, end = const

        #we no longer fold the whole sequence but only the constraint region +- window size
        if window is not None:
            if dist:
                tostart = start - (window - abs(dist))
                toend = fend + (window - abs(dist)) + 1
                log.debug(logid+'Considering distance '+str(dist)+' between constraints '+str(const)+' for window extraction from '+str(tostart)+' to '+str(toend))
            else:
                tostart = start - window
                toend = end + window + 1

            if tostart < 0:
                tostart = 0
            if toend > len(fa.seq):
                toend = len(fa.seq)
            if fend is not None and toend < fend:
                log.warning(logid+'Constraint '+str(cons)+' out of sequence range '+str(toend)+'!Skipping!')  # One of the constraints is outside the sequence window
                return
            seqtofold = str(fa.seq[tostart:toend])

        log.debug(logid+','.join([str(start),str(end),str(tostart),str(toend),'SEQUENCE: '+seqtofold]))

        if len(seqtofold) < window:
            log.warning(logid+'Sequence of '+goi+' to short, seqlenght '+str(len(seqtofold))+' with window size '+str(window)+'!Skipping! ')
            raise Exception('Sequence of '+goi+' to short, seqlenght '+str(len(seqtofold))+' with window size '+str(window))

        cstart = start-tostart
        cend =  end-tostart

        if cstart < 0 or cend > len(fa.seq):
            log.warning(logid+'start of constraint '+str(cstart)+' end of constraint '+str(cend)+' while length of sequence '+str(len(fa.seq))+'! Skipping!')
            return

        checklist = []

        checklist.append((cstart,cend))
        if fstart and fend:
            cfstart = fstart-tostart
            cfend = fend-tostart
            checklist.append((cfstart,cfend))
            checklist.append((cstart,cend,cfstart,cfend))

        #data
        data = {'seq' : seqtofold, 'stru' : []}

        #set model details
        md = RNA.md()
        md.max_bp_span = span

        log.debug(logid+'Constraints for '+goi+' are '+str(checklist))

        for check in checklist:
            s,e,os,oe = [None,None,None,None]
            if genecoords:
                if goi in genecoords:
                    gs, ge, gstrand = get_location(genecoords[goi][0])
                    if gstrand != strand:
                        log.warning(logid+'Strand values differ between Gene annotation and FASTA file! Please check your input for '+str(goi))
                else:
                    gs, ge, gstrand = 0, 0, '.'
                    log.warning(logid+'No coords found for gene '+goi+'! Assuming coordinates are already local!')
            else:
                gs = ge = 0
                gstrand = '.'
                log.warning(logid+'No coords found for gene '+goi+'! Assuming coordinates are already local!')

            log.debug(logid+'GENECOORDS: '+str(goi)+': '+str.join(',',[str(gs),str(ge),str(gstrand)]))

            if len(check) < 3:
                s,e = check
                sp, ep = [s+tostart+gs-1, e+tostart+gs]
                gtostart, gtoend = [tostart+gs, toend+gs]
                printcons = str.join('|',[str.join('-',[str(tostart), str(toend)]), str.join('-',[str(gtostart), str(gtoend)]), str.join('-',[str(s), str(e)]), str.join('-',[str(sp), str(ep)])])
            else:
                s,e,os,oe = check
                sp,ep,osp,oep = [s+tostart+gs-1, e+tostart+gs, os+tostart+gs-1, oe+tostart+gs]
                gtostart, gtoend = [tostart+gs, toend+gs]
                log.debug(logid+'PAIRED:'+';'.join(map(str,[s,e,os,oe,gs,ge,sp,ep,osp,oep,tostart,toend,gtostart,gtoend])))
                printcons = str.join('|',[str.join('-', [str(tostart), str(toend)]), str.join('-',[str(gtostart), str(gtoend)]), str.join(':',[str.join('-',[str(s),str(e)]), str.join('-',[str(os),str(oe)])]), str.join(':',[str.join('-',[str(sp),str(ep)]), str.join('-',[str(osp),str(oep)])])])

            # create new fold_compound object
            fc = RNA.fold_compound(data['seq'], md)
            fc_u = RNA.fold_compound(data['seq'], md)
            fc_p = RNA.fold_compound(data['seq'], md)

            # call pf and prop calculation
            gibbs = fc.pf()[1]
            bppm = get_bppm(fc.bpp(), cstart, cend)

            if bppm is None:
                log.error(logid+'Empty bpp matrix returned, stopping here!')
                #                sys.exit(logid+'Empty bpp matrix returned, stopping here!')
                return

            #enforce paired
            fc_p = constrain_paired(fc_p, s, e)
            #enforce unpaired
            fc_u = constrain_unpaired(fc_u, s, e)
            #calculate probs and nrg
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

            outlist.append([fa, fn, gibbs, '0', nrg, printcons, str(window), str(span), outdir, 'unconstraint'])  # unconstraint
            outlist.append([fa, fn, gibbs_u, dg_u, nrg_u, printcons, str(window), str(span), outdir, 'constraint_unpaired'])  # constraint_paired
            outlist.append([fa, fn, gibbs_p, dg_p, nrg_p, printcons, str(window), str(span), outdir, 'constraint_paired'])  # constraint_unpaired

            if os and oe:
                log.debug(logid+'Second constraint: '+str(','.join(map(str,[goi,data['seq'],len(data['seq']),s,e,os,oe]))))
                #enforce both constraints
                #enforce both paired
                fc_p = constrain_paired(fc_p, os, oe)
                #enforce both unpaired
                fc_u = constrain_unpaired(fc_u, os, oe)
                #calculate probs and nrg
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

                outlist.append([fa, fn, gibbs_u, dg_u, nrg_u, printcons, str(window), str(span), outdir, 'bothconstraint_unpaired'])  # bothconstraint_unpaired
                outlist.append([fa, fn, gibbs_p, dg_p, nrg_p, printcons, str(window), str(span), outdir, 'bothconstraint_paired'])  # bothconstraint_paired

                #enforce second constraint
                # First clear old constraints
                fc_u = RNA.fold_compound(data['seq'], md)
                fc_p = RNA.fold_compound(data['seq'], md)
                #enforce second paired
                fc_p = constrain_paired(fc_p, os, oe)
                #enforce second unpaired
                fc_u = constrain_unpaired(fc_u, os, oe)
                #calculate probs and nrg
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

                outlist.append([fa, fn, gibbs_u, dg_u, nrg_u, printcons, str(window), str(span), outdir, 'secondconstraint_unpaired'])  # secondconstraint_unpaired
                outlist.append([fa, fn, gibbs_p, dg_p, nrg_p, printcons, str(window), str(span), outdir, 'secondconstraint_paired'])  # secondconstraint_paired

        write_out(outlist)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


def constrain_temp(fa, temp, window, span, an, save, outdir, queue=None, configurer=None, level=None):
    # refresh model details
    logid = scriptname+'.constrain_temp: '
    log.info(logid+'Constraining Temp to ' + temp)
    try:
        if queue and level:
            configurer(queue, level)

	    md = RNA.md()
	    md.max_bp_span = span
	    md.window_size = window
	    #set temperature
	    md.temperature = temp
	    #create new fold_compound objects
	    fc_t = RNA.fold_compound(str(fa.seq), md, RNA.OPTION_WINDOW)
	    #data
	    data = {'seq' : str(fa.seq), 'stru' : '', 'nrg' : ''}
	    #set model details
	    md = RNA.md()
	    # create new fold_compound object
	    fc = RNA.fold_compound(data['seq'], md, RNA.OPTION_PF)
	    gibbs = fc.pf()
	    bppm = fc.bpp()

	    data['stru'] = gibbs[0]
	    data['nrg'] = gibbs[1]

	    # call prop window calculation
	    log.debug(logid+';'.join(map(str,[gibbs, constrain, conslength])))

	    for item in bppm:
	        for i in range(int(constrain),int(constrain)+conslength):
	            log.debug(bppm.index(item), i, item[i])

	    at = up_to_array(data_t['up'], None, len(fa.seq))

	    diff_nt = an - at

	    #####Set temp, and rewrite save for temp
	    if save:
	        log.info(logid+'SAVINGTEMP')
	        write_temp(fa, str(temp), data_t, diff_nt, str(window), outdir)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def foldaround(seq, fc, pos, clength, gibbs, nrg, queue=None, configurer=None, level=None):
    # here we take the already constraint fc and constrain regions of length clength around it to see what happens at the original binding site

    logid = scriptname+'.foldaround: '
    try:
        if queue and level:
            configurer(queue, level)

        cstart = pos
        cend = pos+clength-1
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
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


def fold_unconstraint(seq, queue=None, configurer=None, level=None):

    logid = scriptname+'.fold_unconstraint: '
    try:
        if queue and level:
            configurer(queue, level)

        # set model details
        md = RNA.md()
        # create new fold_compound object
        fc = RNA.fold_compound(seq, md, RNA.OPTION_PF)
        # call prop window calculation
        gibbs_uc = fc.pf()[1]

        return(gibbs_uc)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


def write_out(resultlist):

    try:
        for result in resultlist:
            fa, fname, gibbs, ddg, nrg, const, window, span, outdir, condition = result
            logid = scriptname+'.write_out: '
            goi, chrom, strand = idfromfa(fa.id)
            temp_outdir = os.path.join(outdir,goi)

            if fname != 'STDOUT':
                if not os.path.exists(temp_outdir):
                    os.makedirs(temp_outdir)
                if not os.path.exists(os.path.join(temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span])+'.gz')):
                    o = gzip.open(os.path.join(temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span])+'.gz'), 'wb')
                    o.write(bytes(str.join('\t',['Condition','FreeNRG(gibbs)','deltaG','OpeningNRG','Constraint'])+'\n',encoding='UTF-8'))
                    o.write(bytes(str.join('\t',[condition, str(gibbs), str(ddg), str(nrg), str(const)])+'\n',encoding='UTF-8'))
                else:
                    log.info(os.path.join(temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span])+'.gz')+' exists, will append!')
                    o = gzip.open(os.path.join(temp_outdir, '_'.join([goi, chrom, strand, fname, const, window, span])+'.gz'), 'ab')
                    o.write(bytes(str.join('\t',[condition, str(gibbs), str(ddg), str(nrg), str(const)])+'\n',encoding='UTF-8'))
            else:
                print(str.join('\t',['Condition','FreeNRG(gibbs)','deltaG','OpeningNRG','Constraint']))
                print(str.join('\t',[condition, str(gibbs), str(ddg), str(nrg), str(const)]))

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


def main(args):

    logid = scriptname+'.main: '
    try:
        queue, listener, worker_configurer = rissmed_logging_setup(args.logdir, args.loglevel, SCRIPTNAME)

        log.info(logid + 'Running ' + SCRIPTNAME + ' on ' + str(args.procs) + ' cores.')
        log.info(logid+'CLI: '+sys.argv[0]+' '+'{}'.format(' '.join([shlex.quote(s) for s in sys.argv[1:]])))

        run_settings, outdir = preprocess(args.sequence, args.constrain, args.conslength, args.outdir, args.genes)
        fold(args.sequence, args.window, args.span, args.unconstraint, args.unpaired, args.paired, args.length, args.gc, args.number, args.constrain, args.conslength, args.alphabet, args.save, args.procs, args.vrna, args.temprange, args.genecoords, outdir, run_settings, queue=queue, configurer=worker_configurer, level=args.loglevel)
                # def preprocess(queue, configurer, level, sequence, window, span, unconstraint, unpaired, paired, length, gc, number, constrain, conslength, alphabet, save, procs, vrna, temprange, outdir, genes, verbosity=False, pattern=None, cutoff=None):
                #preprocess(queue, worker_configurer, loglevel, args.sequence, args.window, args.span, args.unconstraint, args.unpaired, args.paired, args.length, args.gc, args.number, args.constrain, args.conslength, args.alphabet, args.save, args.procs, args.vrna, args.temprange, args.outdir, args.genes, args.verbosity, args.pattern, args.cutoff)

        queue.put_nowait(None)
        listener.join()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


####################
####    MAIN    ####
####################

if __name__ == '__main__':
    scriptname = os.path.basename(__file__).replace('.py', '')
    logid = scriptname+'.main: '
    try:
        args = parseargs_foldcons()
        main(args)

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

# ConstraintFold.py ends here
