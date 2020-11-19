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
import os.path
import pickle
import stdb
from obspy.clients.fdsn import Client
from obspy.geodetics.base import gps2dist_azimuth as epi
from obspy.geodetics import kilometer2degrees as k2d
from obspy.core import Stream, UTCDateTime
from obstools.atacr import utils, EventStream
from pathlib import Path

from argparse import ArgumentParser
from os.path import exists as exist
from numpy import nan

# Main function

def get_event_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_download_event.py` that accompany this package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <Station Database>",
        description="Script used " +
        "to download and pre-process four-component " +
        "(H1, H2, Z and P), two-hour-long seismograms for " +
        "individual events on which to apply the de-noising " +
        "algorithms. Data are requested from the internet using " +
        "the client services framework for a given date range. " +
        "The stations are processed one by one and the data are " +
        "stored to disk.")
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
        help="Specify a comma separated list of station keys " +
        "for which to perform the analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the "
        "dictionary. For instance, providing IU will match with " +
        "all stations in the IU network [Default processes " +
        "all stations in the database]")
    parser.add_argument(
        "-C", "--channels",
        action="store",
        type=str,
        dest="channels",
        default="",
        help="Specify a comma-separated list of channels for " +
        "which to perform the transfer function analysis. " +
        "Possible options are H (for horizontal channels) or P " +
        "(for pressure channel). Specifying H allows " +
        "for tilt correction. Specifying P allows for compliance " +
        "correction. [Default looks for both horizontal and " +
        "pressure and allows for both tilt AND compliance corrections]")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")

    # Server Settings
    ServerGroup = parser.add_argument_group(
        title="Server Settings",
        description="Settings associated with which "
        "datacenter to log into.")
    ServerGroup.add_argument(
        "-S", "--Server",
        action="store",
        type=str,
        dest="Server",
        default="IRIS",
        help="Specify the server to connect to. Options include: BGR, " +
        "ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU, NCEDC, NEIP, " +
        "NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. [Default IRIS]")
    ServerGroup.add_argument(
        "-U", "--User-Auth",
        action="store",
        type=str,
        dest="UserAuth",
        default="",
        help="Enter your IRIS Authentification Username and Password " +
        "(--User-Auth='username:authpassword') to access and download " +
        "restricted data. [Default no user and password]")

    """
    # # Database Settings
    # DataGroup = parser.add_argument_group(parser, title="Local Data Settings", description="Settings associated with defining " \
    #     "and using a local data base of pre-downloaded day-long SAC files.")
    # DataGroup.add_argument("--local-data", action="store", type=str, dest="localdata", default=None, \
    #     help="Specify a comma separated list of paths containing day-long sac files of data already downloaded. " \
    #     "If data exists for a seismogram is already present on disk, it is selected preferentially over downloading " \
    #     "the data using the Client interface")
    # DataGroup.add_argument("--no-data-zero", action="store_true", dest="ndval", default=False, \
    #     help="Specify to force missing data to be set as zero, rather than default behaviour which sets to nan.")
    """

    # Constants Settings
    FreqGroup = parser.add_argument_group(
        title='Frequency Settings',
        description="Miscellaneous frequency settings")
    FreqGroup.add_argument(
        "--sampling-rate",
        action="store",
        type=float,
        dest="new_sampling_rate",
        default=5.,
        help="Specify new sampling rate (float, in Hz). " +
        "[Default 5.]")
    FreqGroup.add_argument(
        "--pre-filt",
        action="store",
        type=str,
        dest="pre_filt",
        default=None,
        help="Specify four comma-separated corner " +
        "frequencies (float, in Hz) for deconvolution " +
        "pre-filter. [Default 0.001,0.005,45.,50.]")

    # Event Selection Criteria
    EventGroup = parser.add_argument_group(
        title="Event Settings",
        description="Settings associated with refining " +
        "the events to include in matching station " +
        "pairs")
    EventGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the event " +
        "search. This will override any station start " +
        "times. [Default start date of each station in " +
        "database]")
    EventGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string " +
        "representing the start time for the event " +
        "search. This will override any station end times " +
        "[Default end date of each station in database]")
    EventGroup.add_argument(
        "--reverse-order", "-R",
        action="store_true",
        dest="reverse",
        default=False,
        help="Reverse order of events. Default behaviour " +
        "starts at oldest event and works towards most " +
        "recent. Specify reverse order and instead the " +
        "program will start with the most recent events " +
        "and work towards older")
    EventGroup.add_argument(
        "--min-mag",
        action="store",
        type=float,
        dest="minmag",
        default=5.5,
        help="Specify the minimum magnitude of event " +
        "for which to search. [Default 5.5]")
    EventGroup.add_argument(
        "--max-mag",
        action="store",
        type=float,
        dest="maxmag",
        default=None,
        help="Specify the maximum magnitude of event " +
        "for which to search. " +
        "[Default None, i.e. no limit]")

    # Geometry Settings
    GeomGroup = parser.add_argument_group(
        title="Geometry Settings",
        description="Settings associatd with the " +
        "event-station geometries")
    GeomGroup.add_argument(
        "--min-dist",
        action="store",
        type=float,
        dest="mindist",
        default=30.,
        help="Specify the minimum great circle distance " +
        "(degrees) between the station and event. " +
        "[Default 30]")
    GeomGroup.add_argument(
        "--max-dist",
        action="store",
        type=float,
        dest="maxdist",
        default=120.,
        help="Specify the maximum great circle distance " +
        "(degrees) between the station and event. " +
        "[Default 120]")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # create channel list
    if len(args.channels) > 0:
        args.channels = args.channels.split(',')
    else:
        args.channels = ['H', 'P']

    for cha in args.channels:
        if cha not in ['H', 'P']:
            parser.error("Error: Channel not recognized " + str(cha))

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                args.endT)
    else:
        args.endT = None

    # Parse User Authentification
    if not len(args.UserAuth) == 0:
        tt = args.UserAuth.split(':')
        if not len(tt) == 2:
            parser.error(
                "Error: Incorrect Username and Password Strings for User " +
                "Authentification")
        else:
            args.UserAuth = tt
    else:
        args.UserAuth = []

    # # Parse Local Data directories
    # if args.localdata is not None:
    #     args.localdata = args.localdata.split(',')
    # else:
    #     args.localdata = []

    # # Check NoData Value
    # if args.ndval:
    #     args.ndval = 0.0
    # else:
    #     args.ndval = nan

    if args.pre_filt is None:
        args.pre_filt = [0.001, 0.005, 45., 50.]
    else:
        args.pre_filt = [float(val) for val in args.pre_filt.split(',')]
        args.pre_filt = sorted(args.pre_filt)
        if (len(args.pre_filt)) != 4:
            raise(Exception(
                "Error: --pre-filt should contain 4 " +
                "comma-separated floats"))

    return args


def main(args=None):

    if args is None:
        # Run Input Parser
        args = get_daylong_arguments()

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

        # Define path to see if it exists
        eventpath = Path('EVENTS') / Path(stkey)
        if not eventpath.is_dir():
            print('Path to '+str(eventpath)+' doesn`t exist - creating it')
            eventpath.mkdir(parents=True)

        # Establish client
        if len(args.UserAuth) == 0:
            client = Client(args.Server)
        else:
            client = Client(
                args.Server, user=args.UserAuth[0], password=args.UserAuth[1])

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
        print("| Searching Possible events:                    |")
        print("|   Start: {0:19s}                  |".format(
            tstart.strftime("%Y-%m-%d %H:%M:%S")))
        print("|   End:   {0:19s}                  |".format(
            tend.strftime("%Y-%m-%d %H:%M:%S")))
        if args.maxmag is None:
            print("|   Mag:   >{0:3.1f}                                 " +
                  "|".format(
                      args.minmag))
        else:
            print(
                "|   Mag:   {0:3.1f} - {1:3.1f}                         " +
                "   |".format(args.minmag, args.maxmag))
        print("| ...                                           |")

        # Get catalogue using deployment start and end
        cat = client.get_events(
            starttime=tstart, endtime=tend,
            minmagnitude=args.minmag, maxmagnitude=args.maxmag)

        # Total number of events in Catalogue
        nevK = 0
        nevtT = len(cat)
        print(
            "|  Found {0:5d} possible events                  " +
            "|".format(nevtT))
        ievs = range(0, nevtT)

        # Select order of processing
        if args.reverse:
            ievs = range(0, nevtT)
        else:
            ievs = range(nevtT-1, -1, -1)

        # Read through catalogue
        for iev in ievs:

            # Extract event
            ev = cat[iev]

            window = 7200.
            new_sampling_rate = 5.

            time = ev.origins[0].time
            dep = ev.origins[0].depth
            lon = ev.origins[0].longitude
            lat = ev.origins[0].latitude
            epi_dist, az, baz = epi(lat, lon, sta.latitude, sta.longitude)
            epi_dist /= 1000.
            gac = k2d(epi_dist)
            mag = ev.magnitudes[0].mag
            if mag is None:
                mag = -9.

            # If distance between 85 and 120 deg:
            if (gac > args.mindist and gac < args.maxdist):

                # Display Event Info
                nevK = nevK + 1
                if args.reverse:
                    inum = iev + 1
                else:
                    inum = nevtT - iev + 1
                print(" ")
                print(
                    "****************************************************")
                print(
                    "* #{0:d} ({1:d}/{2:d}):  {3:13s}".format(
                        nevK, inum, nevtT, time.strftime("%Y%m%d_%H%M%S")))
                print(
                    "*   Origin Time: " + time.strftime("%Y-%m-%d %H:%M:%S"))
                print(
                    "*   Lat: {0:6.2f}; Lon: {1:7.2f}".format(lat, lon))
                print(
                    "*   Dep: {0:6.2f}; Mag: {1:3.1f}".format(dep/1000., mag))
                print(
                    "*     {0:5s} -> Ev: {1:7.2f} km; {2:7.2f} deg; " +
                    "{3:6.2f}; {4:6.2f}".format(
                        sta.station, epi_dist, gac, baz, az))

                t1 = time
                t2 = t1 + window

                # Time stamp
                tstamp = str(time.year).zfill(4)+'.' + \
                    str(time.julday).zfill(3)+'.'
                tstamp = tstamp + str(time.hour).zfill(2) + \
                    '.'+str(time.minute).zfill(2)

                # Define file names (to check if files already exist)
                filename = eventpath / (tstamp+'.event.pkl')
                file1 = eventpath / (tstamp+'.1.SAC')
                file2 = eventpath / (tstamp+'.2.SAC')
                fileZ = eventpath / (tstamp+'.Z.SAC')
                fileP = eventpath / (tstamp+'.P.SAC')

                print()
                print("* Channels selected: " +
                      str(args.channels)+' and vertical')

                # If data file exists, continue
                if filename.exists():
                    if not args.ovr:
                        print("*")
                        print("*   "+str(filename))
                        print("*   -> File already exists, continuing")
                        continue

                if "P" not in args.channels:

                    # Number of channels
                    ncomp = 3

                    # Comma-separated list of channels for Client
                    channels = sta.channel.upper() + '1,' + \
                        sta.channel.upper() + '2,' + \
                        sta.channel.upper() + 'Z'

                    # Get waveforms from client
                    try:
                        print("*   "+tstamp +
                              "                                     ")
                        print("*   -> Downloading Seismic data... ")
                        sth = client.get_waveforms(
                            network=sta.network, station=sta.station,
                            location=sta.location[0], channel=channels,
                            starttime=t1, endtime=t2, attach_response=True)
                        print("*      ...done")
                    except:
                        print(
                            " Error: Unable to download ?H? components - " +
                            "continuing")
                        continue

                    # Make sure length is ok
                    llZ = len(sth.select(component='Z')[0].data)
                    ll1 = len(sth.select(component='1')[0].data)
                    ll2 = len(sth.select(component='2')[0].data)

                    if (llZ != ll1) or (llZ != ll2):
                        print(" Error: lengths not all the same - continuing")
                        continue

                    ll = int(window*sth[0].stats.sampling_rate)

                    if np.abs(llZ - ll) > 1:
                        print(" Error: Time series too short - continuing")
                        print(np.abs(llZ - ll))
                        continue

                elif "H" not in args.channels:

                    # Number of channels
                    ncomp = 2

                    # Comma-separated list of channels for Client
                    channels = sta.channel.upper() + 'Z'

                    # Get waveforms from client
                    try:
                        print("*   "+tstamp +
                              "                                     ")
                        print("*   -> Downloading Seismic data... ")
                        sth = client.get_waveforms(
                            network=sta.network, station=sta.station,
                            location=sta.location[0], channel=channels,
                            starttime=t1, endtime=t2, attach_response=True)
                        print("*      ...done")
                    except:
                        print(
                            " Error: Unable to download ?H? components - " +
                            "continuing")
                        continue
                    try:
                        print("*   -> Downloading Pressure data...")
                        stp = client.get_waveforms(
                            network=sta.network, station=sta.station,
                            location=sta.location[0], channel='??H',
                            starttime=t1, endtime=t2, attach_response=True)
                        print("*      ...done")
                    except:
                        print(" Error: Unable to download ??H component - " +
                              "continuing")
                        continue

                    # Make sure length is ok
                    llZ = len(sth.select(component='Z')[0].data)
                    llP = len(stp[0].data)

                    if (llZ != llP):
                        print(" Error: lengths not all the same - continuing")
                        continue

                    ll = int(window*stp[0].stats.sampling_rate)

                    if np.abs(llZ - ll) > 1:
                        print(" Error: Time series too short - continuing")
                        print(np.abs(llZ - ll))
                        continue

                else:

                    # Comma-separated list of channels for Client
                    ncomp = 4

                    # Comma-separated list of channels for Client
                    channels = sta.channel.upper() + '1,' + \
                        sta.channel.upper() + '2,' + \
                        sta.channel.upper() + 'Z'

                    # Get waveforms from client
                    try:
                        print("*   "+tstamp +
                              "                                     ")
                        print("*   -> Downloading Seismic data... ")
                        sth = client.get_waveforms(
                            network=sta.network, station=sta.station,
                            location=sta.location[0], channel=channels,
                            starttime=t1, endtime=t2, attach_response=True)
                        print("*      ...done")
                    except:
                        print(
                            " Error: Unable to download ?H? components - " +
                            "continuing")
                        continue
                    try:
                        print("*   -> Downloading Pressure data...")
                        stp = client.get_waveforms(
                            network=sta.network, station=sta.station,
                            location=sta.location[0], channel='??H',
                            starttime=t1, endtime=t2, attach_response=True)
                        print("     ...done")
                    except:
                        print(" Error: Unable to download ??H component - " +
                              "continuing")
                        continue

                    # Make sure length is ok
                    llZ = len(sth.select(component='Z')[0].data)
                    ll1 = len(sth.select(component='1')[0].data)
                    ll2 = len(sth.select(component='2')[0].data)
                    llP = len(stp[0].data)

                    if (llZ != ll1) or (llZ != ll2) or (llZ != llP):
                        print(" Error: lengths not all the same - continuing")
                        continue

                    ll = int(window*sth[0].stats.sampling_rate)

                    if np.abs(llZ - ll) > 1:
                        print(" Error: Time series too short - continuing")
                        print(np.abs(llZ - ll))
                        continue

                # Detrend, filter - seismic data
                sth.detrend('demean')
                sth.detrend('linear')
                sth.filter('lowpass', freq=0.5*args.new_sampling_rate,
                           corners=2, zerophase=True)
                sth.resample(args.new_sampling_rate)

                # Extract traces - Z
                trZ = sth.select(component='Z')[0]
                trZ = utils.update_stats(
                    trZ, sta.latitude, sta.longitude, sta.elevation, 'Z')
                trZ.write(str(fileZ), format='SAC')

                # Extract traces - H
                if "H" in args.channels:
                    tr1 = sth.select(component='1')[0]
                    tr2 = sth.select(component='2')[0]
                    tr1 = utils.update_stats(
                        tr1, sta.latitude, sta.longitude, sta.elevation, '1')
                    tr2 = utils.update_stats(
                        tr2, sta.latitude, sta.longitude, sta.elevation, '2')
                    tr1.write(str(file1), format='SAC')
                    tr2.write(str(file2), format='SAC')

                # Remove responses
                print("*   -> Removing responses - Seismic data")
                sth.remove_response(pre_filt=args.pre_filt, output='DISP')

                if "P" in args.channels:
                    # Detrend, filter - pressure data
                    stp.detrend('demean')
                    stp.detrend('linear')
                    stp.filter('lowpass', freq=0.5 *
                               args.new_sampling_rate, corners=2,
                               zerophase=True)
                    stp.resample(args.new_sampling_rate)

                    # Remove responses
                    print("*   -> Removing responses - Pressure data")
                    stp.remove_response(pre_filt=args.pre_filt)

                    trP = stp[0]
                    trP = utils.update_stats(
                        trP, sta.latitude, sta.longitude, sta.elevation, 'P')
                    trP.write(str(fileP), format='SAC')

                else:
                    stp = Stream()

                # Write out SAC data

                eventstream = EventStream(
                    sta, sth, stp, tstamp, lat, lon, time, window,
                    args.new_sampling_rate, ncomp)

                eventstream.save(filename)


if __name__ == "__main__":

    # Run main program
    main()
