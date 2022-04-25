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
from obstools.atacr import utils, DayNoise
from pathlib import Path

from argparse import ArgumentParser
from os.path import exists as exist
from obspy import UTCDateTime
from numpy import nan


def get_dailyspec_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_daily_spectra.py` that accompany this
    package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <indb>",
        description="Script used "
        "to extract shorter windows from the day-long " +
        "seismograms, calculate the power-spectral " +
        "properties, flag windows for outlier PSDs and " +
        "calculate daily averages of the corresponding " +
        "Fourier transforms. The stations are processed " +
        "one by one and the data are stored to disk. The " +
        "program will look for data saved in the previous " +
        "steps and use all available components.")
    parser.add_argument(
        "indb",
        help="Station Database to process from",
        type=str)

    # General Settings
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys " +
        "for which to perform the analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the " +
        "dictionary. For instance, providing IU will match " +
        "with all stations in the IU network. " +
        "[Default processes all stations in the database]")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Event Selection Criteria
    DaysGroup = parser.add_argument_group(
        title="Time Search Settings",
        description="Time settings associated with " +
        "searching for day-long seismograms")
    DaysGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start day for the data search. " +
        "This will override any station start times. " +
        "[Default start date of each station in database]")
    DaysGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the data search. " +
        "This will override any station end times. " +
        "[Default end date of each station n database]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values " +
        "and settings")
    ConstGroup.add_argument(
        "--window",
        action="store",
        type=float,
        dest="window",
        default=7200.,
        help="Specify window length in seconds. " +
        "Default value is highly recommended. "
        "Program may not be stable for large deviations " +
        "from default value. [Default 7200. (or 2 hours)]")
    ConstGroup.add_argument(
        "--overlap",
        action="store",
        type=float,
        dest="overlap",
        default=0.3,
        help="Specify fraction of overlap between windows. " +
        "[Default 0.3 (or 30 percent)]")
    ConstGroup.add_argument(
        "--minwin",
        action="store",
        type=int,
        dest="minwin",
        default=10,
        help="Specify minimum number of 'good' windows " +
        "in any given day to continue with analysis. " +
        "[Default 10]")
    ConstGroup.add_argument(
        "--freq-band",
        action="store",
        type=str,
        dest="pd",
        default=None,
        help="Specify comma-separated frequency limits " +
        "(float, in Hz) over which to calculate spectral " +
        "features used in flagging the bad windows. " +
        "[Default 0.004,2.0]")
    ConstGroup.add_argument(
        "--tolerance",
        action="store",
        type=float,
        dest="tol",
        default=2.0,
        help="Specify parameter for tolerance threshold. " +
        "If spectrum > std*tol, window is flagged as bad. " +
        "[Default 2.0]")
    ConstGroup.add_argument(
        "--alpha",
        action="store",
        type=float,
        dest="alpha",
        default=0.05,
        help="Specify confidence level for f-test, " +
        "for iterative flagging of windows. " +
        "[Default 0.05, or 95 percent confidence]")
    ConstGroup.add_argument(
        "--raw",
        action="store_true",
        dest="raw",
        default=False,
        help="Raw spectra will be used in calculating " +
        "spectral features for flagging. " +
        "[Default uses smoothed spectra]")
    ConstGroup.add_argument(
        "--no-rotation",
        action="store_false",
        dest="calc_rotation",
        default=True,
        help="Do not rotate horizontal components " +
        "to tilt direction. [Default calculates rotation]")

    # Constants Settings
    FigureGroup = parser.add_argument_group(
        title='Figure Settings',
        description="Flags for plotting figures")
    FigureGroup.add_argument(
        "--figQC",
        action="store_true",
        dest="fig_QC",
        default=False,
        help="Plot Quality-Control figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--debug",
        action="store_true",
        dest="debug",
        default=False,
        help="Plot intermediate steps for debugging. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figAverage",
        action="store_true",
        dest="fig_average",
        default=False,
        help="Plot daily average figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figCoh",
        action="store_true",
        dest="fig_coh_ph",
        default=False,
        help="Plot Coherence and Phase figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--save-fig",
        action="store_true",
        dest="saveplot",
        default=False,
        help="Set this option if you wish to save the figure(s). [Default " +
        "does not save figure]")
    FigureGroup.add_argument(
        "--format",
        action="store",
        type=str,
        dest="form",
        default="png",
        help="Specify format of figure. Can be any one of the valid" +
        "matplotlib formats: 'png', 'jpg', 'eps', 'pdf'. [Default 'png']")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    if args.raw:
        args.smooth = False
    else:
        args.smooth = True

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except Exception:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except Exception:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                args.endT)
    else:
        args.endT = None

    # Check input frequency band
    if args.pd is None:
        args.pd = [0.004, 2.0]
    else:
        args.pd = [float(val) for val in args.pd.split(',')]
        args.pd = sorted(args.pd)
        if (len(args.pd)) != 2:
            raise(Exception(
                "Error: --freq-band should contain 2 " +
                "comma-separated floats"))

    return args


def main(args=None):

    if args is None:
        # Run Input Parser
        args = get_dailyspec_arguments()

    # Load Database
    # stdb>0.1.3
    try:
        db, stkeys = stdb.io.load_db(fname=args.indb, keys=args.stkeys)

    # stdb=0.1.3
    except Exception:
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
        if not datapath.is_dir():
            print("\nPath to "+str(datapath)+" doesn`t exist - continuing")
            continue

        # Path where spectra will be saved
        specpath = Path('SPECTRA') / stkey
        if not specpath.is_dir():
            print("\nPath to "+str(specpath)+" doesn`t exist - creating it")
            specpath.mkdir(parents=True)

        # Path where plots will be saved
        if args.saveplot:
            plotpath = specpath / 'PLOTS'
            if not plotpath.is_dir():
                plotpath.mkdir(parents=True)
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
        print("\n|===============================================|")
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

            print("\n"+"*"*60)
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
