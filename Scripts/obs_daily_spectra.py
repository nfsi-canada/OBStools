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
Station selection is specified by a network and 
station code. The data base is provided in stations_db.pkl as a 
`StDb` dictionary.

Usage
-----

.. code-block::

    $ obs_daily_spectra.py -h
    Usage: obs_daily_spectra.py [options] <station database>

"""

# Import modules and functions
import os
import numpy as np
import pickle
import stdb
from obstools import DayNoise
from obstools import utils, options

def main():

    # Run Input Parser
    (opts, indb) = options.get_dailyspec_options()
    print(opts)

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
            print('Path to '+specpath+' doesn`t exist - creating it')
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

        # Get all components
        trN1, trN2, trNZ, trNP = utils.get_data(datapath, tstart, tend)

        # Window size 
        window = 7200.
        overlap = opts.overlap #0.3

        # minimum numer of windows
        minwin = opts.minwin #10

        # NFFT
        stats = trN1[0].stats

        # Time axis
        taxis = np.arange(0., window, trN1[0].stats.delta)

        # Cycle through available data
        for tr1, tr2, trZ, trP in zip(trN1, trN2, trNZ, trNP):

            year = str(tr1.stats.starttime.year).zfill(4)
            jday = str(tr1.stats.starttime.julday).zfill(3)

            print('Calculating noise spectra for key '+stkey+' and day '+year+'.'+jday)
            tstamp = year+'.'+jday+'.'
            filename = specpath + tstamp + 'spectra.pkl'

            if os.path.exists(filename):
                print('file '+filename+' exists - continuing')
                # continue

            # Initialize instance of DayNoise
            daynoise = DayNoise(tr1, tr2, trZ, trP, window, overlap, key=stkey)

            # Quality control to identify outliers
            daynoise.QC_daily_spectra(fig_QC=opts.fig_QC, debug=opts.debug)

            # Check if we have enough good windows
            nwin = np.sum(daynoise.goodwins)
            if nwin < minwin:
                print('Too few good data segments to calculate day spectra')
                # continue
            else:
                print('{0} good windows. Proceeding...'.format(nwin))

            # Average spectra for good windows
            daynoise.average_daily_spectra(fig_average=opts.fig_average, fig_coh_ph=opts.fig_coh_ph)

            # Save to file
            daynoise.save(filename)


if __name__ == "__main__":

    # Run main program
    main()
