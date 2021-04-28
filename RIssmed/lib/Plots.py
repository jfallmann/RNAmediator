# Plots.py ---
#
# Filename: Plots.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Fri Aug 21 10:22:53 2020 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Fri Aug 21 10:23:22 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 2
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

### IMPORTS
import os
import sys
import inspect
import traceback as tb
import logging

try:
    log = logging.getLogger(__name__)  # use module name
    scriptn = os.path.basename(inspect.stack()[-1].filename).replace('.py', '')
    log.debug('LOGGING IN Plots'+str(scriptn)+str(log)+str(log.handlers))
except Exception:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    print(''.join(tbe.format()),file=sys.stderr)


def plot_data(fa, raw, consu, consp, const, xs, cons, saveas, outdir):
    logid = scriptn+'.plot_data: '
    anime = []
    #define xs for constraint line
    consl = []
    try:
        for x in const:
            consl.append(1.25)
        width = 16/100*len(fa.seq)
        height = 9
        fig = plt.figure(figsize=(width,height),dpi=80)
        ax1 = fig.add_subplot(111)

        ax2 = ax1.twiny()
        #   line, = ax.plot([], [], lw=2)
        plt.title("Blue-- = Unconstraint, Green-. = Unpaired, Red = Paired, Gray = Constraint",y=1.075)
        ax1.set_ylabel('Prob unpaired')
        ax1.set_xlabel('Nucleotides')
        #   plt.xticks(range(0,len(fa.seq)+1),(' '+fa.seq),size='small')
        #add lines to plot
        ax1.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-')
        ax1.set_xlim(0,len(fa.seq)+1)
        ax2.set_xlim(ax1.get_xlim())
        ax1.set_xticks(range(0,len(fa.seq)+1))
        ax1.set_xticklabels((' '+fa.seq), ha="right")
        #   ax2.set_xlabel(r"Modified x-axis: $1/(1+X)$")
        ax2.set_xticks(range(1,len(fa.seq)+1))
        ax2.set_xticklabels(range(1,len(fa.seq)+1), rotation=45, ha="right")
        # We change the fontsize of minor ticks label
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax1.tick_params(axis='both', which='minor', labelsize=4)
        ax2.tick_params(axis='both', which='major', labelsize=5)
        ax2.tick_params(axis='both', which='minor', labelsize=3)
        goi, chrom = fa.id.split(':')[::2]
        strand = str(fa.id.split(':')[3].split('(')[1][0])
        fig.savefig('StruCons_'+goi+'_'+cons+'.'+saveas)
        plt.close()
        #   anime.append(plt.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-'))
        #   return anime
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def plot_temp(fa, raw, temp, xs, saveas, outdir):
    logid = scriptn+'.plot_temp: '
    try:
        anime = []
        #define xs for constraint line
        width = 16/100*len(fa.seq)
        height = 9
        fig = plt.figure(figsize=(width,height),dpi=80)
        ax1 = fig.add_subplot(111)
        plt.title("Blue-- = "+temp+" degree",y=1.075)
        ax1.set_ylabel('Prob unpaired')
        ax1.set_xlabel('Nucleotides')
        #add lines to plot
        ax1.plot(xs, raw, 'b-')
        ax1.set_xlim(0,len(fa.seq)+1)
        ax1.set_xticks(range(0,len(fa.seq)+1))
        ax1.set_xticklabels((' '+fa.seq))
        # We change the fontsize of minor ticks label
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax1.tick_params(axis='both', which='minor', labelsize=4)
        goi, chrom = fa.id.split(':')[::2]
        strand = str(fa.id.split(':')[3].split('(')[1][0])
        fig.savefig('TempCons_'+goi+'_'+temp+'.'+saveas)
        plt.close()
        #   anime.append(plt.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-'))
        #   return anime
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


#
# Plots.py ends here
