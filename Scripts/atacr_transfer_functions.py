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
import sys
import numpy as np
from obspy import UTCDateTime
import pickle
import stdb
from obstools.atacr.classes import StaNoise, Power, Cross, Rotation, TFNoise
from obstools.atacr import utils, plot, options


def main():

    # Run Input Parser
    (opts, indb) = options.get_transfer_options()

    # Load Database
    db = stdb.io.load_db(fname=indb)

    # Construct station key loop
    allkeys = db.keys()
    sorted(allkeys)

    # Extract key subset
    if len(opts.stkeys) > 0:
        stkeys = []
        for skey in opts.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = db.keys()
        sorted(stkeys)

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        if not opts.skip_daily:
            # Path where spectra are located
            specpath = 'SPECTRA/' + stkey + '/'
            if not os.path.isdir(specpath):
                raise(Exception(
                    "Path to "+specpath+" doesn't exist - aborting"))

        if not opts.skip_clean:
            # Path where average spectra will be saved
            avstpath = 'AVG_STA/' + stkey + '/'
            if not os.path.isdir(avstpath):
                print("Path to "+avstpath +
                      " doesn't exist - skipping cleaned station spectra")
                opts.skip_clean = True

        if opts.skip_daily and opts.skip_clean:
            print("skipping both daily and clean spectra")
            continue

        # Path where transfer functions will be located
        tfpath = 'TF_STA/' + stkey + '/'
        if not os.path.isdir(tfpath):
            print("Path to "+tfpath+" doesn't exist - creating it")
            os.makedirs(tfpath)

        # Get catalogue search start time
        if opts.startT is None:
            tstart = sta.startdate
        else:
            tstart = opts.startT

        # Get catalogue search end time
        if opts.endT is None:
            tend = sta.enddate
        else:
            tend = opts.endT

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

        # Filename for output transfer functions
        dstart = str(tstart.year).zfill(4)+'.'+str(tstart.julday).zfill(3)+'-'
        dend = str(tend.year).zfill(4)+'.'+str(tend.julday).zfill(3)+'.'
        fileavst = avstpath + dstart + dend + 'avg_sta.pkl'

        # Find all files in directories
        spectra_files = os.listdir(specpath)
        if not opts.skip_clean:
            average_files = os.listdir(avstpath)

        if not opts.skip_daily:

            day_transfer_functions = []

            # Cycle through available files
            for filespec in spectra_files:

                year = filespec.split('.')[0]
                jday = filespec.split('.')[1]

                print()
                print(
                    "*********************************************" +
                    "***************")
                print("* Calculating transfer functions for key " +
                      stkey+" and day "+year+"."+jday)
                tstamp = year+'.'+jday+'.'
                filename = tfpath + tstamp + 'transfunc.pkl'

                # Load file
                file = open(specpath+filespec, 'rb')
                daynoise = pickle.load(file)
                file.close()

                # Load spectra into TFNoise object
                daytransfer = TFNoise(
                    daynoise.f, daynoise.power,
                    daynoise.cross, daynoise.rotation,
                    daynoise.tf_list)

                # Calculate the transfer functions
                daytransfer.transfer_func()

                # Store the frequency axis
                f = daytransfer.f

                # Append to list of transfer functions
                day_transfer_functions.append(daytransfer.transfunc)

                # Save daily transfer functions to file
                daytransfer.save(filename)

        if not opts.skip_clean:

            # Cycle through available files
            for fileavst in average_files:

                name = fileavst.split('avg_sta')

                print()
                print(
                    "*********************************************" +
                    "***************")
                print("* Calculating transfer functions for key " +
                      stkey+" and range "+name[0])
                filename = tfpath + name[0] + 'transfunc.pkl'

                # Load file
                file = open(avstpath+fileavst, 'rb')
                stanoise = pickle.load(file)
                file.close()

                # Load spectra into TFNoise object - no Rotation object
                # for station averages
                rotation = Rotation(None, None, None)
                statransfer = TFNoise(
                    stanoise.f, stanoise.power,
                    stanoise.cross, rotation, stanoise.tf_list)

                # Calculate the transfer functions
                statransfer.transfer_func()

                # Store the frequency axis
                f = statransfer.f

                # Extract the transfer functions
                sta_transfer_functions = statransfer.transfunc

                # Save average transfer functions to file
                statransfer.save(filename)

        if opts.fig_TF:
            plot.fig_TF(f, day_transfer_functions, daynoise.tf_list,
                        sta_transfer_functions, stanoise.tf_list, skey=stkey,
                        save=opts.saveplot, form=opts.form)


if __name__ == "__main__":

    # Run main program
    main()
