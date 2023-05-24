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
from obspy import UTCDateTime, Stream
import pickle
import stdb
from obstools.atacr import EventStream
from obstools.atacr import utils, plotting
from pathlib import Path

from argparse import ArgumentParser
from os.path import exists as exist
from numpy import nan


def get_correct_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_correct_event.py` that accompanies this
    package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <indb>",
        description="Script used "
        "to extract transfer functions between various " +
        "components, and use them to clean vertical " +
        "component of OBS data for selected events. The " +
        "noise data can be those obtained from the daily " +
        "spectra (i.e., from `obs_daily_spectra.py`) "
        "or those obtained from the averaged noise spectra " +
        "(i.e., from `obs_clean_spectra.py`). Flags are " +
        "available to specify the source of data to use as " +
        "well as the time range for given events. "
        "The stations are processed one by one and the " +
        "data are stored to disk in a new 'CORRECTED' folder.")
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)

    # General Settings
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station " +
        "keys for which to perform the analysis. These must be "
        "contained within the station database. Partial keys " +
        "will be used to match against those in the "
        "dictionary. For instance, providing IU will match with " +
        "all stations in the IU network. [Default processes "
        "all stations in the database]")
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
        "searching for specific event-related seismograms")
    DaysGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start day for the event search. "
        "This will override any station start times. " +
        "[Default start date of each station in database]")
    DaysGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the event search. "
        "This will override any station end times. [Default " +
        "end date of each station in database]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default " +
        "values and settings")
    ConstGroup.add_argument(
        "--skip-daily",
        action="store_true",
        dest="skip_daily",
        default=False,
        help="Skip daily spectral averages in application " +
        "of transfer functions. [Default False]")
    ConstGroup.add_argument(
        "--skip-clean",
        action="store_true",
        dest="skip_clean",
        default=False,
        help="Skip cleaned spectral averages in " +
        "application of transfer functions. " +
        "[Default False]")
    ConstGroup.add_argument(
        "--fmin",
        action="store",
        type=float,
        dest="fmin",
        default="0.006666666666666667",
        help="Low frequency corner (in Hz) for " +
        "plotting the raw (un-corrected) seismograms. "
        "Filter is a 2nd order, zero phase butterworth " +
        "filter. [Default 1./150.]")
    ConstGroup.add_argument(
        "--fmax",
        action="store",
        type=float,
        dest="fmax",
        default="0.1",
        help="High frequency corner (in Hz) for " +
        "plotting the raw (un-corrected) seismograms. "
        "Filter is a 2nd order, zero phase butterworth " +
        "filter. [Default 1./10.]")

    # Constants Settings
    FigureGroup = parser.add_argument_group(
        title='Figure Settings',
        description="Flags for plotting figures")
    FigureGroup.add_argument(
        "--figRaw",
        action="store_true",
        dest="fig_event_raw",
        default=False,
        help="Plot raw seismogram figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figClean",
        action="store_true",
        dest="fig_plot_corrected",
        default=False,
        help="Plot cleaned vertical seismogram figure. " +
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

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except Exception:
            parser.error(
                "Error: Cannot construct UTCDateTime from " +
                "start time: " + args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except Exception:
            parser.error(
                "Error: Cannot construct UTCDateTime from " +
                "end time: " + args.endT)
    else:
        args.endT = None

    if args.skip_clean and args.skip_daily:
        parser.error(
            "Error: cannot skip both daily and clean averages")

    return args


def main(args=None):

    if args is None:
        # Run Input Parser
        args = get_correct_arguments()

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

        # Path where transfer functions will be located
        transpath = Path('TF_STA') / stkey
        if not transpath.is_dir():
            raise(Exception("Path to "+str(transpath) +
                            " doesn`t exist - aborting"))

        # Path where event data are located
        eventpath = Path('EVENTS') / stkey
        if not eventpath.is_dir():
            raise(Exception("Path to "+str(eventpath) +
                            " doesn`t exist - aborting"))

        # Path where plots will be saved
        if args.saveplot:
            plotpath = eventpath / 'PLOTS'
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

        # Get all components
        trE1, trE2, trEZ, trEP = utils.get_event(eventpath, tstart, tend)

        # Find all TF files in directory
        p = list(transpath.glob('*.*'))
        trans_files = [x for x in p if x.is_file()]

        # Check if folders contain anything
        if not trans_files:
            raise(Exception("There are no transfer functions in folder " +
                            str(transpath)))

        # Cycle through available data
        for tr1, tr2, trZ, trP in zip(trE1, trE2, trEZ, trEP):

            eventstream = EventStream(tr1, tr2, trZ, trP)

            # Check if Trace is from SAC file with event info
            evlo = None
            evla = None
            if hasattr(trZ.stats, 'sac'):
                if hasattr(trZ.stats.sac, 'evlo'):
                    evlo = trZ.stats.sac.evlo
                    evla = trZ.stats.sac.evla

            if args.fig_event_raw:
                fname = stkey + '.' + eventstream.tstamp + '.raw'
                plot = plotting.fig_event_raw(
                    eventstream,
                    fmin=args.fmin, fmax=args.fmax)

                if plotpath:
                    plot.savefig(
                        plotpath / (fname + '.' + args.form),
                        dpi=300, bbox_inches='tight', format=args.form)
                else:
                    plot.show()

            # Cycle through corresponding TF files
            for transfile in trans_files:

                # Skip hidden files and folders
                if transfile.name[0] == '.':
                    continue

                tfprefix = transfile.name.split('transfunc')[0]
                print(tfprefix)

                # This case refers to the "cleaned" spectral averages
                if len(tfprefix) > 9:
                    if not args.skip_clean:

                        yr1 = tfprefix.split('-')[0].split('.')[0]
                        jd1 = tfprefix.split('-')[0].split('.')[1]
                        yr2 = tfprefix.split('-')[1].split('.')[0]
                        jd2 = tfprefix.split('-')[1].split('.')[1]
                        date1 = UTCDateTime(yr1+'-'+jd1)
                        date2 = UTCDateTime(yr2+'-'+jd2)
                        dateev = eventstream.evtime

                        if dateev >= date1 and dateev <= date2:
                            print(str(transfile) +
                                  " file found - applying transfer functions")

                            try:
                                file = open(transfile, 'rb')
                                tfaverage = pickle.load(file)
                                file.close()
                            except Exception:
                                print("File "+str(transfile) +
                                      " exists but cannot be loaded")
                                continue

                            # List of possible transfer functions for station
                            # average files
                            eventstream.correct_data(tfaverage)

                            correct_sta = eventstream.correct
                            if args.fig_plot_corrected:
                                fname = eventstream.prefix + '.sta_corrected'
                                plot = plotting.fig_event_corrected(
                                    eventstream, tfaverage.tf_list)
                                # Save or show figure
                                if plotpath:
                                    plot.savefig(
                                        plotpath / (fname + '.' + args.form),
                                        dpi=300, bbox_inches='tight',
                                        format=args.form)
                                else:
                                    plot.show()

                            # Save corrected data to disk
                            correctpath = eventpath / 'CORRECTED'
                            if not correctpath.is_dir():
                                correctpath.mkdir(parents=True)
                            file = correctpath / eventstream.prefix
                            eventstream.save(str(file) + '.sta.pkl')

                            # Now save as SAC files
                            for key, value in tfaverage.tf_list.items():
                                if value and eventstream.ev_list[key]:

                                    # Postfix
                                    nameZ = '.sta.' + key + '.'
                                    nameZ += sta.channel + 'Z.SAC'

                                    # Add Prefix and Postfix
                                    fileZ = str(file) + nameZ

                                    # Select Z component and update trace
                                    trZ = eventstream.trZ.copy()
                                    trZ.data = eventstream.correct[key]
                                    trZ = utils.update_stats(
                                        trZ, sta.latitude, sta.longitude,
                                        sta.elevation, sta.channel+'Z',
                                        evla=evla,
                                        evlo=evlo)

                                    # Save as SAC file
                                    trZ.write(str(fileZ), format='SAC')

                # This case refers to the "daily" spectral averages
                else:
                    if not args.skip_daily:
                        evprefix = eventstream.tstamp.split('.')
                        evstamp = evprefix[0]+'.'+evprefix[1]+'.'
                        if tfprefix == evstamp:
                            print(str(transfile) +
                                  " file found - applying transfer functions")

                            try:
                                file = open(transfile, 'rb')
                                tfaverage = pickle.load(file)
                                file.close()
                            except Exception:
                                print("File "+str(transfile) +
                                      " exists but cannot be loaded")
                                continue

                            # List of possible transfer functions for station
                            # average files
                            eventstream.correct_data(tfaverage)

                            correct_day = eventstream.correct
                            if args.fig_plot_corrected:
                                fname = eventstream.prefix + '.day_corrected'
                                plot = plotting.fig_event_corrected(
                                    eventstream, tfaverage.tf_list)
                                # Save or show figure
                                if plotpath:
                                    plot.savefig(
                                        plotpath / (fname + '.' + args.form),
                                        dpi=300, bbox_inches='tight',
                                        format=args.form)
                                else:
                                    plot.show()

                            # Save corrected data to disk
                            correctpath = eventpath / 'CORRECTED'
                            if not correctpath.is_dir():
                                correctpath.mkdir(parents=True)
                            file = correctpath / eventstream.prefix
                            eventstream.save(str(file) + '.day.pkl')

                            # Now save as SAC files
                            for key, value in tfaverage.tf_list.items():
                                if value and eventstream.ev_list[key]:

                                    # Postfix
                                    nameZ = '.day.' + key + '.'
                                    nameZ += sta.channel+'Z.SAC'

                                    # Add Prefix and Postfix
                                    fileZ = str(file) + nameZ

                                    # Select Z component and update trace
                                    trZ = eventstream.trZ.copy()
                                    trZ.data = eventstream.correct[key]
                                    trZ = utils.update_stats(
                                        trZ, sta.latitude, sta.longitude,
                                        sta.elevation, sta.channel+'Z',
                                        evla=evla,
                                        evlo=evlo)

                                    # Save as SAC file
                                    trZ.write(str(fileZ), format='SAC')


if __name__ == "__main__":

    # Run main program
    main()
