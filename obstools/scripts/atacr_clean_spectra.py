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
import copy

from obspy import UTCDateTime, read_inventory
from obstools.atacr import StaNoise, Power, Cross, Rotation
from obstools.atacr import utils, plotting

from pathlib import Path
from argparse import ArgumentParser
from os.path import exists as exist


def get_cleanspec_arguments(argv=None):

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
        "--flag-freqs",
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
        "--allFigs",
        action="store_true",
        dest="allfigs",
        default=False,
        help="Plot all figures, except for those created by '--debug'. " +
        "Supercedes all other figure creation arguments. [Default False]")
    FigureGroup.add_argument(
        "--debug",
        action="store_true",
        dest="debug",
        default=False,
        help="Plot intermediate steps for debugging. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figQC",
        action="store_true",
        dest="fig_QC",
        default=False,
        help="Plot Quality-Control figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figAverage",
        action="store_true",
        dest="fig_average",
        default=False,
        help="Plot daily average figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figCross",
        action="store_true",
        dest="fig_av_cross",
        default=False,
        help="Plot cross-spectra figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figTilt",
        action="store_true",
        dest="fig_tilt",
        default=False,
        help="Plot coherence, phase and tilt orientation figure. " +
        "[Default does not plot figure]")
    FigureGroup.add_argument(
        "--figTiltPolar",
        action="store_true",
        dest="fig_tilt_polar",
        default=False,
        help="Plot tilt orientation figure in a polar projection. " +
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

    # Check Extension
    ext = args.indb.split('.')[-1]

    if ext not in ['pkl', 'xml', 'csv']:
        parser.error("Must supply a station list in .pkl, .csv or .xml format ")

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
            raise Exception(
                "Error: --freq-band should contain 2 " +
                "comma-separated floats")

    if args.allfigs:
        args.fig_QC = True
        args.fig_average = True
        args.fig_av_cross = True
        args.fig_tilt = True
        args.fig_tilt_polar = True

    return args


def main(args=None):

    print()
    print("###################################################################")
    print("#       _                                          _              #")
    print("#   ___| | ___  __ _ _ __      ___ _ __   ___  ___| |_ _ __ __ _  #")
    print("#  / __| |/ _ \/ _` | '_ \    / __| '_ \ / _ \/ __| __| '__/ _` | #")
    print("# | (__| |  __/ (_| | | | |   \__ \ |_) |  __/ (__| |_| | | (_| | #")
    print("#  \___|_|\___|\__,_|_| |_|___|___/ .__/ \___|\___|\__|_|  \__,_| #")
    print("#                        |_____|  |_|                             #")
    print("#                                                                 #")
    print("###################################################################")
    print()

    if args is None:
        # Run Input Parser
        args = get_cleanspec_arguments()

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
            raise Exception(
                "Path to " + str(specpath) +
                " doesn`t exist - aborting")

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
        tlocs = copy.copy(sta.location)
        if len(tlocs) == 0:
            tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0:
                tlocs.append("--")

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
        filetilt = avstpath / (dstart+dend+'tilt.csv')

        if fileavst.exists():
            if not args.ovr:
                print("*   -> file "+str(fileavst)+" exists - continuing")
                continue

        # Containers for power and cross spectra
        coh_all = []
        ad_all = []
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

        # Date + tilt list
        date_list = []
        tiltdir_list = []
        tiltang_list = []
        coh_list = []

        # Loop through each day withing time range
        while t1 < tend:

            year = str(t1.year).zfill(4)
            jday = str(t1.julday).zfill(3)

            tstamp = year+'.'+jday
            filespec = specpath / (tstamp + '.spectra.pkl')

            # Load file if it exists
            if filespec.exists():
                print("\n"+"*"*60)
                print('* Calculating noise spectra for key ' +
                      stkey+' and day '+year+'.'+jday)
                print("*   -> file "+str(filespec)+" found - loading")
                file = open(filespec, 'rb')
                daynoise = pickle.load(file)
                file.close()
                tiltdir_list.append(daynoise.rotation.tilt_dir)
                tiltang_list.append(daynoise.rotation.tilt_ang)
                date_list.append(t1.date)
                coh_list.append(daynoise.rotation.coh_value)
                stanoise += daynoise
            else:
                t1 += 3600.*24.
                continue

            coh_all.append(daynoise.rotation.coh)
            ad_all.append(daynoise.rotation.ad)
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
        ad_all = np.array(ad_all)
        ph_all = np.array(ph_all)
        tiltdir_list = np.array(tiltdir_list)
        tiltang_list = np.array(tiltang_list)
        date_list = np.array(date_list)
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
        stanoise.QC_sta_spectra(
            pd=args.pd,
            tol=args.tol,
            alpha=args.alpha,
            fig_QC=args.fig_QC,
            debug=args.debug,
            save=plotpath,
            form=args.form)

        # Average spectra for good days
        stanoise.average_sta_spectra(
            fig_average=args.fig_average,
            save=plotpath,
            form=args.form)

        if args.fig_av_cross:
            fname = stkey + '.' + dstart + dend + 'av_coherence'
            plot = plotting.fig_av_cross(
                stanoise.f,
                coh, stanoise.gooddays,
                'Coherence',
                stanoise.ncomp,
                key=stkey,
                lw=0.5)

            if plotpath:
                plot.savefig(
                    str(plotpath / (fname + '.' + args.form)),
                    dpi=300, bbox_inches='tight', format=args.form)
                plot.close()
            else:
                plot.show()

            fname = stkey + '.' + dstart + dend + 'av_admittance'
            plot = plotting.fig_av_cross(
                stanoise.f,
                ad,
                stanoise.gooddays,
                'Admittance',
                stanoise.ncomp,
                key=stkey,
                lw=0.5)

            if plotpath:
                plot.savefig(
                    str(plotpath / (fname + '.' + args.form)),
                    dpi=300, bbox_inches='tight', format=args.form)
                plot.close()
            else:
                plot.show()

            fname = stkey + '.' + dstart + dend + 'av_phase'
            plot = plotting.fig_av_cross(
                stanoise.f,
                ph,
                stanoise.gooddays,
                'Phase',
                stanoise.ncomp,
                key=stkey,
                marker=',',
                lw=0)

            if plotpath:
                plot.savefig(
                    str(plotpath / (fname + '.' + args.form)),
                    dpi=300, bbox_inches='tight', format=args.form)
                plot.close()
            else:
                plot.show()

        if args.fig_tilt and stanoise.phi is not None:
            fname = stkey + '.' + dstart + dend + 'tilt_date'
            plot = plotting.fig_tilt_date(
                stanoise.gooddays,
                coh_all,
                ph_all,
                ad_all,
                stanoise.phi,
                tiltdir_list,
                tiltang_list,
                date_list)

            if plotpath:
                plot.savefig(
                    str(plotpath / (fname + '.' + args.form)),
                    dpi=300, bbox_inches='tight', format=args.form)
                plot.close()
            else:
                plot.show()

        if args.fig_tilt_polar and stanoise.phi is not None:
            fname = stkey + '.' + dstart + dend + 'tilt_polar'
            plot = plotting.fig_tilt_polar_date(
                stanoise.gooddays,
                coh_all,
                ph_all,
                ad_all,
                stanoise.phi,
                tiltdir_list,
                tiltang_list,
                date_list)

            if plotpath:
                plot.savefig(
                    str(plotpath / (fname + '.' + args.form)),
                    dpi=300, bbox_inches='tight', format=args.form)
                plot.close()
            else:
                plot.show()

        # Save to file
        print()
        print("* Clean station spectra saved to: ")
        print("*   "+str(fileavst))
        stanoise.save(fileavst)

        # Write out events
        print()
        print("* Tilt orientation and coherence as function of time saved to: ")
        print("*   "+str(filetilt))
        print()
        fid = open(filetilt, 'w')
        fid.writelines("Date, Tilt dir. (deg CW from H1), " +
                       "Tilt ang. (deg CW from H1), Max coherence\n")
        for i in range(len(tiltdir_list)):
            line1 = "{0},{1:4.1f},{2:4.2f},{3:4.2f}\n".format(
                date_list[i], tiltdir_list[i], tiltang_list[i], coh_list[i])
            fid.writelines(line1.replace(" ", ""))
        fid.close()


if __name__ == "__main__":

    # Run main program
    main()
