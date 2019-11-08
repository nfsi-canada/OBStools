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
Program obs_download_data.py
----------------------------

Downloads four-component (H1, H2, Z and P), day-long seismograms 
to use in noise corrections of vertical
component data. Station selection is specified by a network and 
station code. The data base is provided in stations_db.pkl as a 
`StDb` dictionary.

Usage
-----

.. code-block::

    $ obs_download_data.py -h
    Usage: obs_download_data.py [options] <station database>

    Script used to download and pre-process four-component (H1, H2, Z and P), 
    day-long seismograms to use in noise corrections of vertical component data. 
    This version requests data on the fly for a given date
    range. Data are requested from the internet using the client services framework.
    The stations are processed one by one and the data are stored to disk. 

    Options:
      -h, --help            show this help message and exit
      --keys=STKEYS         Specify a comma separated list of station keys for
                            which to perform the analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -v, -V, --verbose     Specify to increase verbosity.
      -O, --overwrite       Force the overwriting of pre-existing data.
                            Default behaviour prompts for those that already
                            exist. Selecting overwrite and skip (ie, both flags)
                            negate each other, and both are set to false (every
                            repeat is prompted). [Default False]
      -K, --skip-existing   Skip any event for which existing data
                            exist on disk. Selecting skip and overwrite (ie, both flags)
                            negate each other, and both are set to False (every
                            repeat is prompted). [Default False]

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

      Local Data Settings:
        Settings associated with defining and using a local data base of pre-
        downloaded day-long SAC files.

        --local-data=LOCALDATA
                            Specify a comma separated list of paths containing
                            day-long sac files of data already downloaded. If data
                            exists for a seismogram is already present on disk, it
                            is selected preferentially over downloading the data
                            using the Client interface
        --no-data-zero      Specify to force missing data to be set as zero,
                            rather than default behaviour which sets to nan.

      Event Settings:
        Settings associated with refining the events to include in matching
        station pairs

        --start-time=STARTT
                            Specify a UTCDateTime compatible string representing
                            the start time for the event search. This will
                            override any station start times. [Default more recent
                            start date for each station pair]
        --end-time=ENDT     Specify a UTCDateTime compatible string representing
                            the start time for the event search. This will
                            override any station end times [Default older end date
                            for each the pair of stations]
        -R, --reverse-order
                            Reverse order of events. Default behaviour starts at
                            oldest event and works towards most recent. Specify
                            reverse order and instead the program will start with
                            the most recent events and work towards older
        --min-mag=MINMAG    Specify the minimum magnitude of event for which to
                            search. [Default 6.0]
        --max-mag=MAXMAG    Specify the maximum magnitude of event for which to
                            search. [Default None, i.e. no limit]

      Geometry Settings:
        Settings associatd with the event-station geometries

        --min-dist=MINDIST  Specify the minimum great circle distance (degrees)
                            between the station and event. [Default 85]
        --max-dist=MAXDIST  Specify the maximum great circle distance (degrees)
                            between the station and event. [Default 120]

      Parameter Settings:
        Miscellaneous default values and settings

        --Vp=VP             Specify default P velocity value. [Default 6.0 km/s]
        --SNR=MSNR          Specify the SNR threshold used to determine whether
                            events are processedc. [Default 7.5]
        --window=DTS        Specify time window length before and after the SKS
                            arrival. The total window length is 2*dst. [Default
                            120 s]
        --max-delay=MAXDT   Specify the maximum delay time. [Default 4 s]
        --time-increment=DDT
                            Specify the time increment. [Default 0.1 s]
        --angle-increment=DPHI
                            Specify the angle increment. [Default 1 d]
        --transverse-SNR=SNRTLIM
                            Specify the minimum SNR Threshold for the Transverse
                            component to be considered Non-Null. [Default 1.]
"""

# Import modules and functions
import numpy as np
import os.path
import pickle
import glob
import stdb
from obspy import UTCDateTime
from obspy.core import AttribDict
from obspy.clients.fdsn import Client
from obstools import utils, options

# Main function
def main():

    # Run Input Parser
    (opts, indb) = options.get_daylong_options()
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

        # Define path to see if it exists
        datapath = 'DATA/' + stkey + '/'
        if not os.path.isdir(datapath): 
            print('Path to '+datapath+' doesn`t exist - creating it')
            os.makedirs(datapath)

        # Establish client
        if len(opts.UserAuth) == 0:
            client = Client(opts.Server)
        else:
            client = Client(opts.Server, user=opts.UserAuth[0], password=opts.UserAuth[1])

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
        print("|      Start time: {0:19s}          |".format(sta.startdate.strftime("%Y-%m-%d")))
        print("|      End time:   {0:19s}          |".format(sta.enddate.strftime("%Y-%m-%d")))
        print("|-----------------------------------------------|")
        print("| Searching day-long files:                     |")
        print("|   Start: {0:19s}                  |".format(tstart.strftime("%Y-%m-%d")))
        print("|   End:   {0:19s}                  |".format(tend.strftime("%Y-%m-%d")))
        print("| ...                                           |")


        # Split into 24-hour long segments
        dt = 3600.*24.
        new_sampling_rate = opts.new_sampling_rate

        t1 = tstart
        t2 = tstart + dt

        while t2 <= tend:

            # Time stamp
            tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'
     
            # Define file names (to check if files already exist)
            file1 = datapath + tstamp + '.' + sta.channel + '1.SAC'
            file2 = datapath + tstamp + '.' + sta.channel + '2.SAC'
            fileZ = datapath + tstamp + '.' + sta.channel + 'Z.SAC'
            fileP = datapath + tstamp + '.' + sta.channel + 'H.SAC'

            # If data file exists, continue
            if glob.glob(fileZ) and glob.glob(file1) and glob.glob(file2) and glob.glob(fileP): 
                if not opts.ovr:
                    print("| "+tstamp+"*SAC                                 |")
                    print("| -> Files already exist, continuing            |")
                    t1 += dt
                    t2 += dt
                    continue

            channels = sta.channel.upper() + '1,' + sta.channel.upper() + '2,' + sta.channel.upper() + 'Z'

            # Get waveforms from client
            try:
                print("| "+tstamp+"*SAC                                 |")
                print("| -> Downloading Seismic data... ")
                sth = client.get_waveforms(network=sta.network, station=sta.station, location=sta.location[0], \
                        channel=channels, starttime=t1, endtime=t2, attach_response=True)
                print("     ...done")
            except:
                print(" Error: Unable to download ?H? components - continuing")
                t1 += dt
                t2 += dt
                continue
            try:
                print("| -> Downloading Pressure data...")
                stp = client.get_waveforms(network=sta.network, station=sta.station, location=sta.location[0], \
                        channel='??H', starttime=t1, endtime=t2, attach_response=True)
                print("     ...done")
            except:
                print(" Error: Unable to download ??H component - continuing")
                t1 += dt
                t2 += dt
                continue

            # Make sure length is ok
            llZ = len(sth.select(component='Z')[0].data)
            ll1 = len(sth.select(component='1')[0].data)
            ll2 = len(sth.select(component='2')[0].data)
            llP = len(stp[0].data)

            if (llZ != ll1) or (llZ != ll2) or (llZ != llP):
                print(" Error: lengths not all the same - continuing")
                t1 += dt
                t2 += dt
                continue

            ll = int(dt*sth[0].stats.sampling_rate)

            if np.abs(llZ - ll) > 1:
                print(" Error: Time series too short - continuing")
                print(np.abs(llZ - ll))
                t1 += dt
                t2 += dt
                continue

            # Remove responses
            print("| -> Removing responses - Seismic data")
            sth.remove_response(pre_filt=opts.pre_filt, output='DISP')
            print("| -> Removing responses - Pressure data")
            stp.remove_response(pre_filt=opts.pre_filt)

            # Detrend, filter - seismic data
            sth.detrend('demean')
            sth.detrend('linear')
            sth.filter('lowpass', freq=0.5*new_sampling_rate, corners=2, zerophase=True)
            sth.resample(new_sampling_rate)

            # Detrend, filter - pressure data
            stp.detrend('demean')
            stp.detrend('linear')
            stp.filter('lowpass', freq=0.5*new_sampling_rate, corners=2, zerophase=True)
            stp.resample(new_sampling_rate)

            # Extract traces
            trZ = sth.select(component='Z')[0]
            tr1 = sth.select(component='1')[0]
            tr2 = sth.select(component='2')[0]
            trP = stp[0]

            # Update stats
            tr1 = utils.update_stats(tr1, sta.latitude, sta.longitude, sta.elevation)
            tr2 = utils.update_stats(tr2, sta.latitude, sta.longitude, sta.elevation)
            trZ = utils.update_stats(trZ, sta.latitude, sta.longitude, sta.elevation)
            trP = utils.update_stats(trP, sta.latitude, sta.longitude, sta.elevation)

            # Save traces
            tr1.write(file1, format='sac')
            tr2.write(file2, format='sac')
            trZ.write(fileZ, format='sac')
            trP.write(fileP, format='sac')

            t1 += dt
            t2 += dt


if __name__ == "__main__":

    # Run main program
    main()
