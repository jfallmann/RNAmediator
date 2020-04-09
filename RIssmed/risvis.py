#!/usr/bin/env python3
# risvis.py ---
#
# Filename: risvis.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Apr  7 16:51:00 2020 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Thu Apr  9 15:22:38 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 123
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

from flask import Flask, flash, redirect, render_template, request,   session, abort,send_from_directory,send_file,jsonify
import pandas as pd
import json
import sys,os,gzip
import lib.logger
import inspect
import traceback as tb
import argparse
import shlex

from lib.logger import makelogdir, setup_logger
# Create log dir
makelogdir('logs')
scriptname=os.path.basename(__file__)
logname = scriptname
log = setup_logger(name=logname, log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M')

class DataStore():
    GeneID=None
    Constraint=None
    Zscore=None
    data=None

def parseargs():
    parser = argparse.ArgumentParser(description='Visualize base pairing probabilties before and after constraining.')
    parser.add_argument("-f", "--file", type=str, help='Result file to visualize')
    parser.add_argument("--loglevel", type=str, default='INFO', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def readData(file):
    try:
        logid = scriptname+'.readData: '
        log.info('READING INPUT FILE')
        return pd.read_csv(file, delimiter='\t', names=['Chr','Start','End','Constraint','Accessibility_difference','Strand','Distance_to_constraint','Accessibility_no_constraint','Accessibility_constraint','Energy_Difference','Kd_change','Zscore'])#,'ChrBS','StartBS','EndBS','NameBS','ScoreBS','StrandBS'])

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def serve(file, app):
    try:
        logid = scriptname+'.serve: '
        df = readData(file)
        log.info(df)
        data=DataStore()

        @app.route('/',methods=['GET','POST'])
        def prep_json():
            # Get values from form
            data.GeneID = request.form.get('GeneID', None)
            data.Constraint = request.form.get('Constraint', None)
            data.Zscore = request.form.get('Zscore', '0')

            # Filter
            GeneID = data.GeneID if hasattr(data, 'GeneID') and data.GeneID is not None and data.GeneID != 'GeneID' and data.GeneID != '' else str(df.loc[0,'Constraint']).split('|')[0]
            Constraint = data.Constraint if hasattr(data, 'Constraint') and data.Constraint is not None and data.Constraint != 'GeneID' and data.Constraint != '' else str(df.loc[0,'Constraint']).split('|')[2]
            Zscore = data.Zscore if hasattr(data, 'Zscore') else 0

            log.info(str([GeneID, Constraint, Zscore]))
            plotdata=None

            if GeneID is not None and GeneID != 'GeneID' and GeneID != '':
                # Filter the data frame (df)
                subdf = df[df['Constraint'].str.contains(GeneID) & df['Constraint'].str.contains(Constraint)]
                subdf = df[['Start', 'Accessibility_difference', 'Accessibility_no_constraint', 'Accessibility_constraint', 'Zscore']]

                d = {"name": GeneID, "values": []}

                for line in subdf.values:
                    pos = line[0]
                    diff = line[1]
                    prob_no = line[2]
                    prob_cons = line[3]
                    zscore = line[4]

                    d['values'].append([pos, diff, prob_no, prob_cons, zscore])

                #Dump data to json
                dump = json.dumps(d)#Save to datastore
                with open(GeneID+'.json','w') as f:
                    f.write(dump)
                data.data = json.loads(dump)#Save to temporary variable
                plotdata = data.data

            return render_template("stats.html",GeneID=GeneID,Constraint=Constraint,Zscore=Zscore,data=plotdata)

        @app.route('/get-data',methods=['GET','POST'])
        def returnData():
            log.info('get: '+str(data.data))
            return jsonify(data.data)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        log = setup_logger(name=logname, log_file='logs/'+logname, logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
        log.setLevel(args.loglevel)

        log.info(logid+'Running '+scriptname)
        log.info(logid+'CLI: '+sys.argv[0]+'{}'.format(' '.join( [shlex.quote(s) for s in sys.argv[1:]] )))

        log.info('STARTING FLASK SERVICE')

        app = Flask(__name__, static_url_path='', static_folder='vis/static', template_folder='vis/templates')
        serve(args.file, app)
        app.run(debug=True)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


#
# risvis.py ends here
