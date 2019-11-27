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
Program obs_daily_spectra.py
----------------------------

Extracts two-hour-long windows from the day-long seismograms, calculates 
power-spectral densities and flags windows for outlier from the PSD properties. 
Station selection is specified by a network and station code. The data base 
is provided as a `StDb` dictionary.

Usage
-----

.. code-block::

    $ obs_daily_spectra.py -h
    Usage: obs_daily_spectra.py [options] <station database>

    Script used to extract two-hour-long windows from the day-long seismograms,
    calculate the power-spectral properties, flag windows for outlier PSDs and
    calculate daily averages of the corresponding Fourier transforms. The stations
    are processed one by one and the data are stored to disk.

    Options:
      -h, --help           show this help message and exit
      --keys=STKEYS        Specify a comma separated list of station keys for
                           which to perform the analysis. These must be contained
                           within the station database. Partial keys will be used
                           to match against those in the dictionary. For instance,
                           providing IU will match with all stations in the IU
                           network. [Default processes all stations in the
                           database]
      -O, --overwrite      Force the overwriting of pre-existing data. [Default
                           False]

      Parameter Settings:
        Miscellaneous default values and settings

        --overlap=OVERLAP  Specify fraction of overlap between two-hour-long
                           windows. [Default 0.3 (or 30%)]
        --minwin=MINWIN    Specify minimum number of 'good' windows in any given
                           day to continue with analysis. [Default 10]
        --freq-band=PD     Specify comma-separated frequency limits (float, in Hz)
                           over which to calculate spectral features used in
                           flagging the days/windows. [Default 0.004,2.0]
        --tolerance=TOL    Specify parameter for tolerance threshold. If spectrum
                           > std*tol, window is flagged as bad. [Default 1.5]
        --alpha=ALPHA      Specify confidence level for f-test, for iterative
                           flagging of windows. [Default 0.05, or 95% confidence]
        --raw              Raw spectra will be used in calculating spectral
                           features for flagging. [Default uses smoothed spectra]
        --no-rotation      Do not rotate horizontal components to tilt direction.
                           [Default calculates rotation]

      Figure Settings:
        Flags for plotting figures

        --figQC            Plot Quality-Control figure. [Default does not plot
                           figure]
        --debug            Plot intermediate steps for debugging. [Default does
                           not plot figure]
        --figAverage       Plot daily average figure. [Default does not plot
                           figure]
        --figCoh           Plot Coherence and Phase figure. [Default does not plot
                           figure]

      Time Search Settings:
        Time settings associated with searching for day-long seismograms

        --start=STARTT     Specify a UTCDateTime compatible string representing
                           the start day for the data search. This will override
                           any station start times. [Default start date of each
                           station in database]
        --end=ENDT         Specify a UTCDateTime compatible string representing
                           the start time for the event search. This will override
                           any station end times. [Default end date of each
                           station n database]

"""

# Import modules and functions
import os
import numpy as np
import pickle
import stdb
from obstools import DayNoise
from obstools.atacr import utils, options

def main():

    # Run Input Parser
    (opts, indb) = options.get_dailyspec_options()

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

        # Path where data are located
        datapath = 'DATA/' + stkey + '/'

        # Path where spectra will be saved
        specpath = 'SPECTRA/' + stkey + '/'
        if not os.path.isdir(specpath): 
            print("Path to "+specpath+" doesn`t exist - creating it")
            os.makedirs(specpath)

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

        # Get all components
        trN1, trN2, trNZ, trNP = utils.get_data(datapath, tstart, tend)

        # Window size 
        window = opts.window
        overlap = opts.overlap #0.3

        # minimum numer of windows
        minwin = opts.minwin #10

        # NFFT
        stats = trN1[0].stats

        # Time axis
        taxis = np.arange(0., window, trNZ[0].stats.delta)

        # Cycle through available data
        for tr1, tr2, trZ, trP in zip(trN1, trN2, trNZ, trNP):

            year = str(trZ.stats.starttime.year).zfill(4)
            jday = str(trZ.stats.starttime.julday).zfill(3)

            print(" ")
            print("***************************************************************")
            print("* Calculating noise spectra for key "+stkey+" and day "+year+"."+jday)
            tstamp = year+'.'+jday+'.'
            filename = specpath + tstamp + 'spectra.pkl'

            if os.path.exists(filename):
                if not opts.ovr:
                    print("   -> file "+filename+" exists - continuing")
                    continue

            # Initialize instance of DayNoise
            daynoise = DayNoise(tr1, tr2, trZ, trP, window, overlap, key=stkey)

            # Quality control to identify outliers
            daynoise.QC_daily_spectra(pd=opts.pd, tol=opts.tol, alpha=opts.alpha, smooth=opts.smooth, 
                fig_QC=opts.fig_QC, debug=opts.debug)

            # Check if we have enough good windows
            nwin = np.sum(daynoise.goodwins)
            if nwin < minwin:
                print("*   Too few good data segments to calculate average day spectra")
                # continue
            else:
                print("*   {0} good windows. Proceeding...".format(nwin))

            # Average spectra for good windows
            daynoise.average_daily_spectra(calc_rotation=opts.calc_rotation, fig_average=opts.fig_average, 
                fig_coh_ph=opts.fig_coh_ph)

            # Save to file
            daynoise.save(filename)


if __name__ == "__main__":

    # Run main program
    main()
