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
Program obs_transfer_functions.py
---------------------------------

Calculates transfer functions using the noise windows flagged as "good", for either
a single day (from `obs_daily_spectra.py`) or for those averaged over several days
(from `obs_clean_spectra.py`), if available. The transfer functions are stored to disk.

Station selection is specified by a network and 
station code. The data base is provided as a 
`StDb` dictionary.

Usage
-----

.. code-block::

    $ obs_transfer_functions.py -h
    Usage: obs_transfer_functions.py [options] <station database>

    Script used to calculate transfer functions between various components, to be
    used in cleaning vertical component of OBS data. The noise data can be those
    obtained from the daily spectra (i.e., from `obs_daily_spectra.py`) or those
    obtained from the averaged noise spectra (i.e., from `obs_clean_spectra.py`).
    Flags are available to specify the source of data to use as well as the time
    range over which to calculate the transfer functions. The stations are
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

        --skip-daily    Skip daily spectral averages in construction of transfer
                        functions. [Default False]
        --skip-clean    Skip cleaned spectral averages in construction of transfer
                        functions. Defaults to True if data cannot be found in
                        default directory. [Default False]

      Figure Settings:
        Flags for plotting figures

        --figTF         Plot transfer function figure. [Default does not plot
                        figure]

      Time Search Settings:
        Time settings associated with searching for day-long seismograms

        --start=STARTT  Specify a UTCDateTime compatible string representing the
                        start day for the data search. This will override any
                        station start times. [Default start date of each station
                        in database]
        --end=ENDT      Specify a UTCDateTime compatible string representing the
                        start time for the event search. This will override any
                        station end times. [Default end date of each station in
                        database]

"""

# Import modules and functions
import os
import sys
import numpy as np
from obspy import UTCDateTime
import pickle
from obstools import StaNoise, Power, Cross, Rotation, TFNoise
from obstools import utils, plot, options

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
            if not os.path.isdir(avstpath): 
                raise(Exception("Path to "+specpath+" doesn`t exist - aborting"))

        if not opts.skip_clean:
            # Path where average spectra will be saved
            avstpath = 'AVG_STA/' + stkey + '/'
            if not os.path.isdir(avstpath): 
                print("Path to "+avstpath+" doesn`t exist - skipping averaged spectra over multiple days")
                opts.skip_clean = True

        # Path where transfer functions will be located
        tfpath = 'TF_STA/' + stkey + '/'
        if not os.path.isdir(tfpath): 
            print("Path to "+tfpath+" doesn`t exist - creating it")
            os.makedirs(tfpath)

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


        # Filename for output transfer functions
        dstart = str(tstart.year).zfill(4)+'.'+str(tstart.julday).zfill(3)+'-'
        dend = str(tend.year).zfill(4)+'.'+str(tend.julday).zfill(3)+'.'
        fileavst = avstpath + dstart + dend + 'avg_sta.pkl'


        # Find all files in directories
        spectra_files = os.listdir(specpath)
        if not opts.skip_clean:
            average_files = os.listdir(avstpath)

        if not opts.skip_daily:

            # List of possible transfer functions for Daily files
            TF_list_day = {'ZP': True, 'Z1':True, 'Z2-1':True, 'ZP-21':True, 'ZH':True, 'ZP-H':True}

            day_transfer_functions = []

            # Cycle through available files
            for filespec in spectra_files:

                year = filespec.split('.')[0]
                jday = filespec.split('.')[1]

                print('Calculating transfer functions for key '+stkey+' and day '+year+'.'+jday)
                tstamp = year+'.'+jday+'.'
                filename = tfpath + tstamp + 'transfunc.pkl'

                # Load file
                file = open(specpath+filespec, 'rb')
                daynoise = pickle.load(file)
                file.close()

                # Load spectra into TFNoise object
                daytransfer = TFNoise(daynoise.f, daynoise.power, daynoise.cross, daynoise.rotation, TF_list_day)

                # Calculate the transfer functions
                daytransfer.transfer_func()

                # Store the frequency axis
                f = daytransfer.f

                # Append to list of transfer functions
                day_transfer_functions.append(daytransfer.transfunc)

                # Save daily transfer functions to file
                daytransfer.save(filename)

        if not opts.skip_clean:

            # List of possible transfer functions for station average files
            TF_list_sta = {'ZP': True, 'Z1':True, 'Z2-1':True, 'ZP-21':True, 'ZH':False, 'ZP-H':False}

            # Cycle through available files
            for fileavst in average_files:

                name = fileavst.split('avg_sta')

                print('Calculating transfer functions for key '+stkey+' and range '+name[0])
                filename = tfpath + name[0] + 'transfunc.pkl'

                # Load file
                file = open(avstpath+fileavst, 'rb')
                stanoise = pickle.load(file)
                file.close()

                # Load spectra into TFNoise object - no Rotation object for station averages
                rotation = Rotation(None, None, None)
                statransfer = TFNoise(stanoise.f, stanoise.power, stanoise.cross, rotation, TF_list_sta)

                # Calculate the transfer functions
                statransfer.transfer_func()

                # Store the frequency axis
                f = statransfer.f

                # Extract the transfer functions
                sta_transfer_functions = statransfer.transfunc

                # Save average transfer functions to file
                statransfer.save(filename)

        if opts.fig_TF:
            plot.fig_TF(f, day_transfer_functions, sta_transfer_functions, key=stkey)

if __name__ == "__main__":

    # Run main program
    main()

