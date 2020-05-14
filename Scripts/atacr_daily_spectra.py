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
import pickle
import stdb
from obstools.atacr import DayNoise
from obstools.atacr import utils, arguments
from pathlib import Path


def main():

    # Run Input Parser
    args = arguments.get_dailyspec_arguments()

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

        # Path where data are located
        datapath = Path('DATA') / stkey

        # Path where spectra will be saved
        specpath = Path('SPECTRA') / stkey
        if not specpath.is_dir():
            print()
            print("Path to "+str(specpath)+" doesn`t exist - creating it")
            specpath.mkdir()

        # Path where plots will be saved
        if args.saveplot:
            plotpath = specpath / 'PLOTS'
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
        print()
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

        # Get all components
        trN1, trN2, trNZ, trNP = utils.get_data(datapath, tstart, tend)

        # Window size
        window = args.window
        overlap = args.overlap

        # minimum numer of windows
        minwin = args.minwin

        # Time axis
        taxis = np.arange(0., window, trNZ[0].stats.delta)

        # Cycle through available data
        for tr1, tr2, trZ, trP in zip(trN1, trN2, trNZ, trNP):

            year = str(trZ.stats.starttime.year).zfill(4)
            jday = str(trZ.stats.starttime.julday).zfill(3)

            print()
            print(
                "************************************************************")
            print("* Calculating noise spectra for key " +
                  stkey+" and day "+year+"."+jday)
            tstamp = year+'.'+jday+'.'
            filename = specpath / (tstamp+'spectra.pkl')

            if filename.exists():
                if not args.ovr:
                    print("*   -> file "+str(filename)+" exists - continuing")
                    continue

            # Initialize instance of DayNoise
            daynoise = DayNoise(tr1, tr2, trZ, trP, window, overlap, key=stkey)

            # Quality control to identify outliers
            daynoise.QC_daily_spectra(
                pd=args.pd, tol=args.tol, alpha=args.alpha,
                smooth=args.smooth, fig_QC=args.fig_QC,
                save=plotpath, form=args.form, debug=args.debug)

            # Check if we have enough good windows
            nwin = np.sum(daynoise.goodwins)
            if nwin < minwin:
                print("*   Too few good data segments to calculate " +
                      "average day spectra")
                # continue
            else:
                print("*   {0} good windows. Proceeding...".format(nwin))

            # Average spectra for good windows
            daynoise.average_daily_spectra(
                calc_rotation=args.calc_rotation,
                fig_average=args.fig_average,
                fig_coh_ph=args.fig_coh_ph,
                save=plotpath, form=args.form)

            # Save to file
            daynoise.save(filename)


if __name__ == "__main__":

    # Run main program
    main()
