#!/usr/bin/env python

# Copyright 2019 Pascal Audet & Helen Janiszewski
#
# This file is part of OBStools.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# Import modules and functions
import os
import numpy as np
from obspy import UTCDateTime
import pickle
import stdb
from obstools.atacr import StaNoise, Power, Cross, Rotation, TFNoise
from obstools.atacr import utils, plot, arguments


def main():

    # Run Input Parser
    args = arguments.get_correct_arguments()

    # Load Database
    db = stdb.io.load_db(fname=args.indb)

    # Construct station key loop
    allkeys = db.keys()
    sorted(allkeys)

    # Extract key subset
    if len(args.stkeys) > 0:
        stkeys = []
        for skey in args.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = db.keys()
        sorted(stkeys)

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Path where transfer functions will be located
        transpath = 'TF_STA/' + stkey + '/'
        if not os.path.isdir(transpath):
            raise(Exception("Path to "+transpath+" doesn`t exist - aborting"))

        # Path where event data are located
        eventpath = 'EVENTS/' + stkey + '/'
        if not os.path.isdir(eventpath):
            raise(Exception("Path to "+eventpath+" doesn`t exist - aborting"))

        # Path where plots will be saved
        if args.saveplot: 
            plotpath = eventpath + 'PLOTS/'
            if not os.path.isdir(plotpath):
                os.makedirs(plotpath)
        else:
            plotpath = False

        # Get catalogue search start time
        if args.startT is None:
            tstart = sta.startdate
        else:
            tstart = args.startT

        # Get catalogue search end time
        if args.endT is None:
            tend = sta.enddate
        else:
            tend = args.endT

        if tstart > sta.enddate or tend < sta.startdate:
            continue

        # Temporary print locations
        tlocs = sta.location
        if len(tlocs) == 0:
            tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0:
                tlocs[il] = "--"
        sta.location = tlocs

        # Update Display
        print(" ")
        print(" ")
        print("|===============================================|")
        print("|===============================================|")
        print("|                   {0:>8s}                    |".format(
            sta.station))
        print("|===============================================|")
        print("|===============================================|")
        print("|  Station: {0:>2s}.{1:5s}                            |".format(
            sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}  |".format(
            sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(
            sta.longitude, sta.latitude))
        print("|      Start time: {0:19s}          |".format(
            sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|      End time:   {0:19s}          |".format(
            sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|-----------------------------------------------|")

        # Find all files in directories
        event_files = os.listdir(eventpath)
        trans_files = os.listdir(transpath)

        # Check if folders contain anything
        if not event_files:
            raise(Exception("There are no events in folder " + eventpath))

        if not trans_files:
            raise(Exception("There are no transfer functions in folder " +
                            transpath))

        # Cycle through available files
        for eventfile in event_files:

            # Skip hidden files and folders
            if eventfile[0] == '.':
                continue

            evprefix = eventfile.split('.')
            evstamp = evprefix[0]+'.'+evprefix[1]+'.'

            evDateTime = UTCDateTime(evprefix[0]+'-'+evprefix[1])
            if evDateTime >= tstart and evDateTime <= tend:

                # Load event file
                try:
                    file = open(eventpath+eventfile, 'rb')
                    eventstream = pickle.load(file)
                    file.close()
                except:
                    print("File "+eventpath+eventfile +
                          " exists but cannot be loaded")
                    continue

            else:
                continue

            if args.fig_event_raw:
                fname = stkey + '.' + evstamp + 'raw'
                plot.fig_event_raw(eventstream, fmin=args.fmin, fmax=args.fmax,
                    save=plotpath, fname=fname, form=args.form)

            # Cycle through corresponding TF files
            for transfile in trans_files:

                # Skip hidden files and folders
                if transfile[0] == '.':
                    continue

                tfprefix = transfile.split('transfunc')[0]

                # This case refers to the "cleaned" spectral averages
                if len(tfprefix) > 9:
                    if not args.skip_clean:
                        yr1 = tfprefix.split('-')[0].split('.')[0]
                        jd1 = tfprefix.split('-')[0].split('.')[1]
                        yr2 = tfprefix.split('-')[1].split('.')[0]
                        jd2 = tfprefix.split('-')[1].split('.')[1]
                        date1 = UTCDateTime(yr1+'-'+jd1)
                        date2 = UTCDateTime(yr2+'-'+jd2)
                        dateev = UTCDateTime(evprefix[0]+'-'+evprefix[1])
                        if dateev >= date1 and dateev <= date2:
                            print(transpath+transfile +
                                  " file found - applying transfer functions")

                            try:
                                file = open(transpath+transfile, 'rb')
                                tfaverage = pickle.load(file)
                                file.close()
                            except:
                                print("File "+transpath+transfile +
                                      " exists but cannot be loaded")
                                continue

                            # List of possible transfer functions for station
                            # average files
                            eventstream.correct_data(tfaverage)

                            correct = eventstream.correct
                            if args.fig_plot_corrected:
                                fname = stkey + '.' + evstamp + 'sta_corrected'
                                plot.fig_event_corrected(
                                    eventstream, tfaverage.tf_list,
                                    save=plotpath, fname=fname, form=args.form)

                # This case refers to the "daily" spectral averages
                else:
                    if not args.skip_daily:
                        if tfprefix == evstamp:
                            print(transpath+transfile +
                                  " file found - applying transfer functions")

                            try:
                                file = open(transpath+transfile, 'rb')
                                tfaverage = pickle.load(file)
                                file.close()
                            except:
                                print("File "+transpath+transfile +
                                      " exists but cannot be loaded")
                                continue

                            # List of possible transfer functions for station
                            # average files
                            eventstream.correct_data(tfaverage)

                            correct = eventstream.correct
                            if args.fig_plot_corrected:
                                fname = stkey + '.' + evstamp + 'day_corrected'
                                plot.fig_event_corrected(
                                    eventstream, tfaverage.tf_list,
                                    save=plotpath, fname=fname, form=args.form)


if __name__ == "__main__":

    # Run main program
    main()
