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
import numpy as np
from obspy import UTCDateTime
import pickle
import stdb
from obstools.atacr import StaNoise, Power, Cross, Rotation, TFNoise
from obstools.atacr import utils, plotting, arguments
from pathlib import Path


def main():

    # Run Input Parser
    args = arguments.get_transfer_arguments()

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

        if not args.skip_daily:
            # Path where spectra are located
            specpath = Path('SPECTRA') / stkey
            if not specpath.is_dir():
                raise(Exception(
                    "Path to "+str(specpath)+" doesn't exist - aborting"))

        if not args.skip_clean:
            # Path where average spectra will be saved
            avstpath = Path('AVG_STA') / stkey
            if not avstpath.is_dir():
                print("Path to "+str(avstpath) +
                      " doesn't exist - skipping cleaned station spectra")
                args.skip_clean = True

        if args.skip_daily and args.skip_clean:
            print("skipping both daily and clean spectra")
            continue

        # Path where transfer functions will be located
        tfpath = Path('TF_STA') / stkey
        if not tfpath.is_dir():
            print("Path to "+str(tfpath)+" doesn't exist - creating it")
            tfpath.mkdir()

        # Path where plots will be saved
        if args.saveplot:
            plotpath = tfpath / 'PLOTS'
            if not plotpath.is_dir():
                plotpath.mkdir()
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

        # Filename for output transfer functions
        dstart = str(tstart.year).zfill(4)+'.'+str(tstart.julday).zfill(3)+'-'
        dend = str(tend.year).zfill(4)+'.'+str(tend.julday).zfill(3)+'.'
        fileavst = avstpath / (dstart+dend+'avg_sta.pkl')

        # Find all files in directories
        p = specpath.glob('*spectra.pkl')
        spectra_files = [x for x in p if x.is_file()]
        if not args.skip_clean:
            p = avstpath.glob('*avg_sta.pkl')
            average_files = [x for x in p if x.is_file()]

        if not args.skip_daily:

            day_transfer_functions = []

            # Cycle through available files
            for filespec in spectra_files:

                year = filespec.name.split('.')[0]
                jday = filespec.name.split('.')[1]

                print()
                print(
                    "*********************************************" +
                    "***************")
                print("* Calculating transfer functions for key " +
                      stkey+" and day "+year+"."+jday)
                tstamp = year+'.'+jday+'.'
                filename = tfpath / (tstamp + 'transfunc.pkl')

                # Load file
                file = open(filespec, 'rb')
                daynoise = pickle.load(file)
                file.close()

                # Load spectra into TFNoise object
                daytransfer = TFNoise(daynoise)

                # Calculate the transfer functions
                daytransfer.transfer_func()

                # Store the frequency axis
                f = daytransfer.f

                # Append to list of transfer functions
                day_transfer_functions.append(daytransfer.transfunc)

                # Save daily transfer functions to file
                daytransfer.save(filename)

        if not args.skip_clean:

            # Cycle through available files
            for fileavst in average_files:

                name = fileavst.name.split('avg_sta')

                print()
                print(
                    "*********************************************" +
                    "***************")
                print("* Calculating transfer functions for key " +
                      stkey+" and range "+name[0])
                filename = tfpath / (name[0] + 'transfunc.pkl')

                # Load file
                file = open(fileavst, 'rb')
                stanoise = pickle.load(file)
                file.close()

                # Load spectra into TFNoise object - no Rotation object
                # for station averages
                rotation = Rotation(None, None, None)
                statransfer = TFNoise(stanoise)

                # Calculate the transfer functions
                statransfer.transfer_func()

                # Store the frequency axis
                f = statransfer.f

                # Extract the transfer functions
                sta_transfer_functions = statransfer.transfunc

                # Save average transfer functions to file
                statransfer.save(filename)

        if args.fig_TF:
            fname = stkey + '.' + 'transfer_functions'
            plot = plotting.fig_TF(f, day_transfer_functions, daynoise.tf_list,
                        sta_transfer_functions, stanoise.tf_list, skey=stkey)

            if plotpath:
                plot.savefig(plotpath / (fname + '.' + args.form),
                            dpi=300, bbox_inches='tight', format=args.form)
            else:
                plot.show()



if __name__ == "__main__":

    # Run main program
    main()
