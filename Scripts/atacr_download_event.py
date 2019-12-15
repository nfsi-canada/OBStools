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

"""
Program atacr_download_event.py
-------------------------------

Downloads up to four-component (H1, H2, Z and P), two-hour-long seismograms 
for individual seismic events to use in noise corrections of vertical
component data. Station selection is specified by a network and 
station code. The data base is provided in stations_db.pkl as a 
`StDb` dictionary.

Usage
-----

.. code-block::

    $ atacr_download_event.py -h
    Usage: atacr_download_event.py [options] <station database>

    Script used to download and pre-process four-component (H1, H2, Z and P), two-
    hour-long seismograms for individual events on which to apply the de-noising
    algorithms. Data are requested from the internet using the client services
    framework for a given date range. The stations are processed one by one and
    the data are stored to disk.

    Options:
      -h, --help            show this help message and exit
      --keys=STKEYS         Specify a comma separated list of station keys for
                            which to perform the analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -C CHANNELS, --channels=CHANNELS
                            Specify a comma-separated list of channels for which
                            to perform the transfer function analysis. Possible
                            options are H (for horizontal channels) or P (for
                            pressure channel). Specifying H allows for tilt
                            correction. Specifying P allows for compliance
                            correction. [Default looks for both horizontal and
                            pressure and allows for both tilt AND compliance
                            corrections]
      -O, --overwrite       Force the overwriting of pre-existing data. [Default
                            False]

      Server Settings:
        Settings associated with which datacenter to log into.

        -S SERVER, --Server=SERVER
                            Specify the server to connect to. Options include:
                            BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU,
                            NCEDC, NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS,
                            USP. [Default IRIS]
        -U USERAUTH, --User-Auth=USERAUTH
                            Enter your IRIS Authentification Username and Password
                            (--User-Auth='username:authpassword') to access and
                            download restricted data. [Default no user and
                            password]

      Event Settings:
        Settings associated with refining the events to include in matching
        station pairs

        --start=STARTT      Specify a UTCDateTime compatible string representing
                            the start time for the event search. This will
                            override any station start times. [Default start date
                            of each station in database]
        --end=ENDT          Specify a UTCDateTime compatible string representing
                            the start time for the event search. This will
                            override any station end times [Default end date of
                            each station in database]
        -R, --reverse-order
                            Reverse order of events. Default behaviour starts at
                            oldest event and works towards most recent. Specify
                            reverse order and instead the program will start with
                            the most recent events and work towards older
        --min-mag=MINMAG    Specify the minimum magnitude of event for which to
                            search. [Default 5.5]
        --max-mag=MAXMAG    Specify the maximum magnitude of event for which to
                            search. [Default None, i.e. no limit]

      Geometry Settings:
        Settings associatd with the event-station geometries

        --min-dist=MINDIST  Specify the minimum great circle distance (degrees)
                            between the station and event. [Default 30]
        --max-dist=MAXDIST  Specify the maximum great circle distance (degrees)
                            between the station and event. [Default 120]

      Frequency Settings:
        Miscellaneous frequency settings

        --sampling-rate=NEW_SAMPLING_RATE
                            Specify new sampling rate (float, in Hz). [Default 5.]
        --pre-filt=PRE_FILT
                            Specify four comma-separated corner frequencies
                            (float, in Hz) for deconvolution pre-filter. [Default
                            0.001,0.005,45.,50.]
                        
"""

# Import modules and functions
import numpy as np
import os.path
import pickle
import glob
import stdb
from obspy.clients.fdsn import Client
from obspy.geodetics.base import gps2dist_azimuth as epi
from obspy.geodetics import kilometer2degrees as k2d
from obspy.core import Stream
from obstools.atacr import utils, options
from obstools.atacr.classes import EventStream

# Main function


def main():

    # Run Input Parser
    (opts, indb) = options.get_event_options()

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

        # Define path to see if it exists
        eventpath = 'EVENTS/' + stkey + '/'
        if not os.path.isdir(eventpath):
            print('Path to '+eventpath+' doesn`t exist - creating it')
            os.makedirs(eventpath)

        # Establish client
        if len(opts.UserAuth) == 0:
            client = Client(opts.Server)
        else:
            client = Client(
                opts.Server, user=opts.UserAuth[0], password=opts.UserAuth[1])

        # Get catalogue search start time
        if opts.startT is None:
            tstart = sta.startdate
        else:
            tstart = opts.startT

        # Get catalogue search end time
        if opts.endT is None:
            tend = sta.enddate
        else:
            tend = opts.endT
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
        if opts.maxmag is None:
            print("|   Mag:   >{0:3.1f}                                 " +
                  "|".format(
                      opts.minmag))
        else:
            print(
                "|   Mag:   {0:3.1f} - {1:3.1f}                         " +
                "   |".format(opts.minmag, opts.maxmag))
        print("| ...                                           |")

        # Get catalogue using deployment start and end
        cat = client.get_events(
            starttime=tstart, endtime=tend,
            minmagnitude=opts.minmag, maxmagnitude=opts.maxmag)

        # Total number of events in Catalogue
        nevK = 0
        nevtT = len(cat)
        print(
            "|  Found {0:5d} possible events                  " +
            "|".format(nevtT))
        ievs = range(0, nevtT)

        # Select order of processing
        if opts.reverse:
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
            if (gac > opts.mindist and gac < opts.maxdist):

                # Display Event Info
                nevK = nevK + 1
                if opts.reverse:
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
                tstamp = str(time.year).zfill(4)+'.' +
                str(time.julday).zfill(3)+'.'
                tstamp = tstamp + str(time.hour).zfill(2) +
                '.'+str(time.minute).zfill(2)

                # Define file names (to check if files already exist)
                filename = eventpath + tstamp + '.event.pkl'

                print()
                print("* Channels selected: " +
                      str(opts.channels)+' and vertical')

                # If data file exists, continue
                if glob.glob(filename):
                    if not opts.ovr:
                        print("*")
                        print("*   "+filename)
                        print("*   -> File already exists, continuing")
                        continue

                if "P" not in opts.channels:

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

                elif "H" not in opts.channels:

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

                # Remove responses
                print("*   -> Removing responses - Seismic data")
                sth.remove_response(pre_filt=opts.pre_filt, output='DISP')
                if "P" in opts.channels:
                    print("*   -> Removing responses - Pressure data")
                    stp.remove_response(pre_filt=opts.pre_filt)

                # Detrend, filter - seismic data
                sth.detrend('demean')
                sth.detrend('linear')
                sth.filter('lowpass', freq=0.5*opts.new_sampling_rate,
                           corners=2, zerophase=True)
                sth.resample(opts.new_sampling_rate)

                if "P" in opts.channels:
                    # Detrend, filter - pressure data
                    stp.detrend('demean')
                    stp.detrend('linear')
                    stp.filter('lowpass', freq=0.5 *
                               opts.new_sampling_rate, corners=2,
                               zerophase=True)
                    stp.resample(opts.new_sampling_rate)
                else:
                    stp = Stream()

                eventstream = EventStream(
                    sta, sth, stp, tstamp, lat, lon, time, window,
                    opts.new_sampling_rate, ncomp)

                eventstream.save(filename)


if __name__ == "__main__":

    # Run main program
    main()
