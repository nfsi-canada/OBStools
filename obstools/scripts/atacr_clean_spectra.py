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
import pickle
import stdb
from obstools.atacr import StaNoise, Power, Cross, Rotation
from obstools.atacr import utils, plotting
from pathlib import Path

from argparse import ArgumentParser
from os.path import exists as exist
from obspy import UTCDateTime
from numpy import nan


def get_cleanspec_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_clean_spectra.py` that accompany this
    package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <indb>",
        description="Script used "
        "to extract daily spectra calculated from " +
        "`obs_daily_spectra.py` and flag days for outlier " +
        "PSDs and calculate spectral averages of the " +
        "corresponding Fourier transforms over the entire " +
        "time period specified. The stations are processed " +
        "one by one and the data are stored to disk.")
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
        "keys for which to perform the analysis. These must " +
        "be contained within the station database. Partial " +
        "keys will be used to match against those in the " +
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
        "[Default end date of each station in database]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values " +
        "and settings")
    ConstGroup.add_argument(
        "--freq-band",
        action="store",
        type=str,
        dest="pd",
        default=None,
        help="Specify comma-separated frequency limits " +
        "(float, in Hz) over which to calculate spectral " +
        "features used in flagging the days/windows. " +
        "[Default 0.004,2.0]")
    ConstGroup.add_argument(
        "--tolerance",
        action="store",
        type=float,
        dest="tol",
        default=1.5,
        help="Specify parameter for tolerance threshold. " +
        "If spectrum > std*tol, window is flagged as bad. " +
        "[Default 1.5]")
    ConstGroup.add_argument(
        "--alpha",
        action="store",
        type=float,
        dest="alpha",
        default=0.05,
        help="Confidence level for f-test, for iterative " +
        "flagging of windows. [Default 0.05, or 95 percent confidence]")

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
        "--figCross",
        action="store_true",
        dest="fig_av_cross",
        default=False,
        help="Plot cross-spectra figure. " +
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
        args = get_cleanspec_arguments()

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

        # Path where spectra are located
        specpath = Path('SPECTRA') / stkey
        if not specpath.is_dir():
            raise(Exception(
                "Path to " + str(specpath) +
                " doesn`t exist - aborting"))

        # Path where average spectra will be saved
        avstpath = Path('AVG_STA') / stkey
        if not avstpath.is_dir():
            print("Path to "+str(avstpath)+" doesn`t exist - creating it")
            avstpath.mkdir(parents=True)

        # Path where plots will be saved
        if args.saveplot:
            plotpath = avstpath / 'PLOTS'
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

        # Filename for output average spectra
        dstart = str(tstart.year).zfill(4)+'.'+str(tstart.julday).zfill(3)+'-'
        dend = str(tend.year).zfill(4)+'.'+str(tend.julday).zfill(3)+'.'
        fileavst = avstpath / (dstart+dend+'avg_sta.pkl')

        if fileavst.exists():
            if not args.ovr:
                print("*   -> file "+str(fileavst)+" exists - continuing")
                continue

        # Containers for power and cross spectra
        coh_all = []
        ph_all = []
        coh_12_all = []
        coh_1Z_all = []
        coh_1P_all = []
        coh_2Z_all = []
        coh_2P_all = []
        coh_ZP_all = []
        ph_12_all = []
        ph_1Z_all = []
        ph_1P_all = []
        ph_2Z_all = []
        ph_2P_all = []
        ph_ZP_all = []
        ad_12_all = []
        ad_1Z_all = []
        ad_1P_all = []
        ad_2Z_all = []
        ad_2P_all = []
        ad_ZP_all = []
        nwins = []

        t1 = tstart

        # Initialize StaNoise object
        stanoise = StaNoise()

        # Loop through each day withing time range
        while t1 < tend:

            year = str(t1.year).zfill(4)
            jday = str(t1.julday).zfill(3)

            tstamp = year+'.'+jday+'.'
            filespec = specpath / (tstamp + 'spectra.pkl')

            # Load file if it exists
            if filespec.exists():
                print("\n"+"*"*60)
                print('* Calculating noise spectra for key ' +
                      stkey+' and day '+year+'.'+jday)
                print("*   -> file "+str(filespec)+" found - loading")
                file = open(filespec, 'rb')
                daynoise = pickle.load(file)
                file.close()
                stanoise += daynoise
            else:
                t1 += 3600.*24.
                continue

            coh_all.append(daynoise.rotation.coh)
            ph_all.append(daynoise.rotation.ph)

            # Coherence
            coh_12_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.c12,
                        daynoise.power.c11,
                        daynoise.power.c22), 50))
            coh_1Z_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.c1Z,
                        daynoise.power.c11,
                        daynoise.power.cZZ), 50))
            coh_1P_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.c1P,
                        daynoise.power.c11,
                        daynoise.power.cPP), 50))
            coh_2Z_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.c2Z,
                        daynoise.power.c22,
                        daynoise.power.cZZ), 50))
            coh_2P_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.c2P,
                        daynoise.power.c22,
                        daynoise.power.cPP), 50))
            coh_ZP_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.cZP,
                        daynoise.power.cZZ,
                        daynoise.power.cPP), 50))

            # Phase
            try:
                ph_12_all.append(
                    180./np.pi*utils.phase(daynoise.cross.c12))
            except Exception:
                ph_12_all.append(None)
            try:
                ph_1Z_all.append(
                    180./np.pi*utils.phase(daynoise.cross.c1Z))
            except Exception:
                ph_1Z_all.append(None)
            try:
                ph_1P_all.append(
                    180./np.pi*utils.phase(daynoise.cross.c1P))
            except Exception:
                ph_1P_all.append(None)
            try:
                ph_2Z_all.append(
                    180./np.pi*utils.phase(daynoise.cross.c2Z))
            except Exception:
                ph_2Z_all.append(None)
            try:
                ph_2P_all.append(
                    180./np.pi*utils.phase(daynoise.cross.c2P))
            except Exception:
                ph_2P_all.append(None)
            try:
                ph_ZP_all.append(
                    180./np.pi*utils.phase(daynoise.cross.cZP))
            except Exception:
                ph_ZP_all.append(None)

            # Admittance
            ad_12_all.append(utils.smooth(utils.admittance(
                daynoise.cross.c12, daynoise.power.c11), 50))
            ad_1Z_all.append(utils.smooth(utils.admittance(
                daynoise.cross.c1Z, daynoise.power.c11), 50))
            ad_1P_all.append(utils.smooth(utils.admittance(
                daynoise.cross.c1P, daynoise.power.c11), 50))
            ad_2Z_all.append(utils.smooth(utils.admittance(
                daynoise.cross.c2Z, daynoise.power.c22), 50))
            ad_2P_all.append(utils.smooth(utils.admittance(
                daynoise.cross.c2P, daynoise.power.c22), 50))
            ad_ZP_all.append(utils.smooth(utils.admittance(
                daynoise.cross.cZP, daynoise.power.cZZ), 50))

            t1 += 3600.*24.

        # Convert to numpy arrays
        coh_all = np.array(coh_all)
        ph_all = np.array(ph_all)
        coh_12_all = np.array(coh_12_all)
        coh_1Z_all = np.array(coh_1Z_all)
        coh_1P_all = np.array(coh_1P_all)
        coh_2Z_all = np.array(coh_2Z_all)
        coh_2P_all = np.array(coh_2P_all)
        coh_ZP_all = np.array(coh_ZP_all)
        ph_12_all = np.array(ph_12_all)
        ph_1Z_all = np.array(ph_1Z_all)
        ph_1P_all = np.array(ph_1P_all)
        ph_2Z_all = np.array(ph_2Z_all)
        ph_2P_all = np.array(ph_2P_all)
        ph_ZP_all = np.array(ph_ZP_all)
        ad_12_all = np.array(ad_12_all)
        ad_1Z_all = np.array(ad_1Z_all)
        ad_1P_all = np.array(ad_1P_all)
        ad_2Z_all = np.array(ad_2Z_all)
        ad_2P_all = np.array(ad_2P_all)
        ad_ZP_all = np.array(ad_ZP_all)

        # Store transfer functions as objects for plotting
        coh = Cross(coh_12_all, coh_1Z_all, coh_1P_all,
                    coh_2Z_all, coh_2P_all, coh_ZP_all)
        ph = Cross(ph_12_all, ph_1Z_all, ph_1P_all,
                   ph_2Z_all, ph_2P_all, ph_ZP_all)
        ad = Cross(ad_12_all, ad_1Z_all, ad_1P_all,
                   ad_2Z_all, ad_2P_all, ad_ZP_all)

        # Quality control to identify outliers
        stanoise.QC_sta_spectra(pd=args.pd, tol=args.tol, alpha=args.alpha,
                                fig_QC=args.fig_QC, debug=args.debug,
                                save=plotpath, form=args.form)

        # Average spectra for good days
        stanoise.average_sta_spectra(
            fig_average=args.fig_average,
            save=plotpath, form=args.form)

        if args.fig_av_cross:
            fname = stkey + '.' + 'av_coherence'
            plot = plotting.fig_av_cross(
                stanoise.f, coh, stanoise.gooddays,
                'Coherence', stanoise.ncomp, key=stkey, lw=0.5)
            # if plotpath.is_dir():
            if plotpath:
                plot.savefig(
                    str(plotpath / (fname + '.' + args.form)),
                    dpi=300, bbox_inches='tight', format=args.form)
            else:
                plot.show()

            fname = stkey + '.' + 'av_admittance'
            plot = plotting.fig_av_cross(
                stanoise.f, ad, stanoise.gooddays,
                'Admittance', stanoise.ncomp, key=stkey, lw=0.5)

            if plotpath:
                plot.savefig(
                    str(plotpath / (fname + '.' + args.form)),
                    dpi=300, bbox_inches='tight', format=args.form)
            else:
                plot.show()

            fname = stkey + '.' + 'av_phase'
            plot = plotting.fig_av_cross(
                stanoise.f, ph, stanoise.gooddays,
                'Phase', stanoise.ncomp, key=stkey, marker=',', lw=0)

            if plotpath:
                plot.savefig(
                    str(plotpath / (fname + '.' + args.form)),
                    dpi=300, bbox_inches='tight', format=args.form)
            else:
                plot.show()

        if args.fig_coh_ph and stanoise.direc is not None:
            fname = stkey + '.' + 'coh_ph'
            plot = plotting.fig_coh_ph(coh_all, ph_all, stanoise.direc)
            if plotpath:
                plot.savefig(
                    str(plotpath / (fname + '.' + args.form)),
                    dpi=300, bbox_inches='tight', format=args.form)
            else:
                plot.show()

        # Save to file
        stanoise.save(fileavst)


if __name__ == "__main__":

    # Run main program
    main()
