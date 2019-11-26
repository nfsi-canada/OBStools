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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Program obs_correct_event.py
----------------------------

Calculates transfer functions using the noise windows flagged as "good", for either
a single day (from `obs_daily_spectra.py`) or for those averaged over several days
(from `obs_clean_spectra.py`), if available. The transfer functions are stored to disk.

Station selection is specified by a network and 
station code. The data base is provided as a 
`StDb` dictionary.

Usage
-----

.. code-block::

    $ obs_correct_event.py -h
    Usage: obs_correct_event.py [options] <station database>

    Script used to extract transfer functions between various components, and use
    them to clean vertical component of OBS data for selected events. The noise
    data can be those obtained from the daily spectra (i.e., from
    `obs_daily_spectra.py`) or those obtained from the averaged noise spectra
    (i.e., from `obs_clean_spectra.py`). Flags are available to specify the source
    of data to use as well as the time range for given events. The stations are
    processed one by one and the data are stored to disk.

    Options:
      -h, --help        show this help message and exit
      --keys=STKEYS     Specify a comma separated list of station keys for which
                        to perform the analysis. These must be contained within
                        the station database. Partial keys will be used to match
                        against those in the dictionary. For instance, providing
                        IU will match with all stations in the IU network.
                        [Default processes all stations in the database]
      -O, --overwrite   Force the overwriting of pre-existing data. [Default
                        False]

      Parameter Settings:
        Miscellaneous default values and settings

        --skip-daily    Skip daily spectral averages in application of transfer
                        functions. [Default False]
        --skip-clean    Skip cleaned spectral averages in application of transfer
                        functions. [Default False]
        --fmin=FMIN     Low frequency corner (in Hz) for plotting the raw (un-
                        corrected) seismograms. Filter is a 2nd order, zero phase
                        butterworth filter. [Default 1./150.]
        --fmax=FMAX     High frequency corner (in Hz) for plotting the raw (un-
                        corrected) seismograms. Filter is a 2nd order, zero phase
                        butterworth filter. [Default 1./10.]

      Figure Settings:
        Flags for plotting figures

        --figRaw        Plot raw seismogram figure. [Default does not plot figure]
        --figClean      Plot cleaned vertical seismogram figure. [Default does not
                        plot figure]

      Time Search Settings:
        Time settings associated with searching for specific event-related
        seismograms

        --start=STARTT  Specify a UTCDateTime compatible string representing the
                        start day for the event search. This will override any
                        station start times. [Default start date of each station
                        in database]
        --end=ENDT      Specify a UTCDateTime compatible string representing the
                        start time for the event search. This will override any
                        station end times. [Default end date of each station in
                        database]

"""

# Import modules and functions
import os
import numpy as np
from obspy import UTCDateTime
import pickle
from obstools import StaNoise, Power, Cross, Rotation, TFNoise
from obstools.atacr import utils, plot, options

def main():

    # Run Input Parser
    (opts, indb) = options.get_correct_options()

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

        # Path where transfer functions will be located
        tfpath = 'TF_STA/' + stkey + '/'
        if not os.path.isdir(tfpath): 
            raise(Exception("Path to "+tfpath+" doesn`t exist - aborting"))

        # Path where event data are located
        eventpath = 'EVENTS/' + stkey + '/'
        if not os.path.isdir(eventpath): 
            raise(Exception("Path to "+tfpath+" doesn`t exist - aborting"))

        # Get catalogue search start time
        if opts.startT is None:
            tstart = sta.startdate
        else:
            tstart = opts.startT

        # Get catalogue search end time
        if opts.endT is None:
            tend = sta.startdate
        else:
            tend = opts.endT

        if tstart > sta.enddate or tend < sta.startdate:
            continue

        # Temporary print locations
        tlocs = sta.location
        if len(tlocs) == 0: tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0: tlocs[il] = "--"
        sta.location = tlocs

        # Update Display
        print(" ")
        print(" ")
        print("|===============================================|")
        print("|===============================================|")
        print("|                   {0:>8s}                    |".format(sta.station))
        print("|===============================================|")
        print("|===============================================|")
        print("|  Station: {0:>2s}.{1:5s}                            |".format(sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}  |".format(sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(sta.longitude, sta.latitude))
        print("|      Start time: {0:19s}          |".format(sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|      End time:   {0:19s}          |".format(sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|-----------------------------------------------|")


        # Find all files in directories
        event_files = os.listdir(eventpath)
        trans_files = os.listdir(transpath)

        # Check if folders contain anything
        if not event_files:
            raise(Exception("There are no events in folder "+eventpath))

        if not trans_files:
            raise(Exception("There are no transfer functions in folder "+transpath))

        # Cycle through available files
        for eventfile in event_files:

            # Skip hidden files and folders
            if eventfile[0]=='.':
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
                    print("File "+eventpath+eventfile+" exists but cannot be loaded")
                    continue

            else:
                continue

            if opts.fig_event_raw:
                plot.fig_event_raw(eventstream, fmin=opts.fmin, fmax=opts.fmax)

            # Cycle through corresponding TF files
            for transfile in trans_files:

                # Skip hidden files and folders
                if transfile[0]=='.':
                    continue

                tfprefix = transfile.split('transfunc')[0]

                # This case refers to the "cleaned" spectral averages
                if len(tfprefix) > 9:
                    if not opts.skip_clean:
                        yr1 = tfprefix.split('-')[0].split('.')[0]
                        jd1 = tfprefix.split('-')[0].split('.')[1]
                        yr2 = tfprefix.split('-')[1].split('.')[0]
                        jd2 = tfprefix.split('-')[1].split('.')[1]
                        if evprefix[0]>=yr1 and evprefix[0] <=yr2:
                            if evprefix[1]>=jd1 and evprefix[1] <=jd2:
                                print(transpath+transfile+" file found - applying transfer functions")

                                try:
                                    file = open(transpath+transfile, 'rb')
                                    tfaverage = pickle.load(file)
                                    file.close()
                                except:
                                    print("File "+transpath+transfile+" exists but cannot be loaded")
                                    continue

                                # List of possible transfer functions for station average files
                                TF_list = tfaverage.tf_list
                                eventstream.correct_data(tfaverage, TF_list)

                                correct = eventstream.correct
                                if opts.plot_corrected:
                                    plot.fig_event_corrected(eventstream, TF_list)

                # This case refers to the "daily" spectral averages
                else:
                    if not opts.skip_daily:
                        if tfprefix==evstamp:
                            print(transpath+transfile+" file found - applying transfer functions")

                            try:
                                file = open(transpath+transfile, 'rb')
                                tfaverage = pickle.load(file)
                                file.close()
                            except:
                                print("File "+transpath+transfile+" exists but cannot be loaded")
                                continue

                            # List of possible transfer functions for station average files
                            TF_list = tfaverage.tf_list
                            eventstream.correct_data(tfaverage, TF_list)

                            correct = eventstream.correct
                            if opts.plot_corrected:
                                plot.fig_event_corrected(eventstream, TF_list)


if __name__ == "__main__":

    # Run main program
    main()

