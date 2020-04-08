# risvis.py ---
#
# Filename: risvis.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Apr  7 16:51:00 2020 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Wed Apr  8 11:01:22 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 19
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
import logging
import inspect
import traceback as tb

def parseargs():
    parser = argparse.ArgumentParser(description='Visualize base pairing probabilties before and after constraining.')
    parser.add_argument("-f", "--file", type=str, help='Result file to visualize')
    parser.add_argument("--loglevel", type=str, default='WARNING', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def serve(file):
    #flask and class definition
    app= Flask(__name__)

    class DataStore():
        GeneID=None
        Constraint=None
        Zscore=None
        prob_no = None
        prob_cons = None
        probdiff = None

    data=DataStore()


    @app.route(“/”,methods=[“GET”,”POST”])

    #Load data
    df = pd.read_csv(file, delimiter='\t', names=['Chr','Start','End','Constraint','Accessibility_difference','Strand','Distance_to_constraint','Accessibility_no_constraint','Accessibility_constraint','Energy_Difference','Kd_change','Zscore','ChrBS','StartBS','EndBS','NameBS','ScoreBS','StrandBS'])

    # Get values from form
    data.GeneID = request.form.get(‘GeneID’,’’)
    data.Constraint = request.form.get(‘Constraint’, '')
    data.Zscore = request.form.get(‘Zscore’, '0')

    def prep_json(df):
        # Filter
        GeneID = data.GeneID if data.GeneId != '' else str.split('|', df.loc[0,'Constraint'])[0]
        Constraint = data.Constraint if data.Constraint != '' else str.split('|', df.loc[0,'Constraint'])[-1]
        Zscore = data.Zscore

        # Filter the data frame (df)
        df = df[GeneID in df.Constraint ]
        df = df[Constraint in df.Constraint]
        df = df[['Start', 'Accessibility_difference', 'Accessibility_no_constraint', 'Accessibility_constraint', 'Zscore']]

        #df1 = df.groupby(['GeneID', 'Start']).sum()
        #df1 = df1.reset_index()#Lets create a dict

        d = {"name": GeneID, "values": []}

        for line in df.values:
            pos = line[0]
            diff = line[1]
            prob_no = line[2]
            prob_cons = line[3]
            zscore = line[4]

            # make a list of keys
            #keys_list = []
            #for item in d['children']:
            #    keys_list.append(item['name'])

            # if 'the_parent' is NOT a key in the flare.json yet, append it
            if not Category in keys_list:
                d['children'].append({"name": Category, "children":      [{"name": Cat, "size": value}]})

        # if 'the_parent' IS a key in the flare.json, add a new child to it
        else:
            d['children'][keys_list.index(Category)]  ['children'].append({"name": Cat, "size": value})

        flare = d

        return render_template(“stats.html”,GeneID=GeneID,Constraint=Constraint,Zscore=Zscore)

    @app.route(“/get-data”,methods=[“GET”,”POST”])

    def returnProdData():
        f=data.Prod
        return jsonify(f)

    app.run(debug=True)

####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        logname = scriptname+'_'+args.constrain
        log = logging.log(name=logname, log_file='logs/'+logname, logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)

        log.info(logid+'Running '+scriptname+' on '+str(args.procs)+' cores.')
        log.info(logid+'CLI: '+sys.argv[0]+'{}'.format(' '.join( [shlex.quote(s) for s in sys.argv[1:]] )))

        serve(args.file)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


#
# risvis.py ends here
