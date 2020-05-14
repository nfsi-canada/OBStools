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
from obstools.atacr import utils, arguments
from pathlib import Path

# Main function


def main():

    # Run Input Parser
    args = arguments.get_daylong_arguments()

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
        datapath = Path('DATA') / Path(stkey)
        if not datapath.is_dir():
            print()
            print('Path to '+str(datapath)+' doesn`t exist - creating it')
            datapath.mkdir()

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
            tend = sta.startdate
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
            sta.startdate.strftime("%Y-%m-%d")))
        print("|      End time:   {0:19s}          |".format(
            sta.enddate.strftime("%Y-%m-%d")))
        print("|-----------------------------------------------|")
        print("| Searching day-long files:                     |")
        print("|   Start: {0:19s}                  |".format(
            tstart.strftime("%Y-%m-%d")))
        print("|   End:   {0:19s}                  |".format(
            tend.strftime("%Y-%m-%d")))

        # Split into 24-hour long segments
        dt = 3600.*24.

        t1 = tstart
        t2 = tstart + dt

        while t2 <= tend:

            # Time stamp
            tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'

            print()
            print(
                "***********************************************************")
            print("* Downloading day-long data for key "+stkey +
                  " and day "+str(t1.year)+"."+str(t1.julday))
            print("*")
            print("* Channels selected: "+str(args.channels)+' and vertical')

            # Define file names (to check if files already exist)
            # Horizontal 1 channel
            file1 = datapath / (tstamp+'.'+sta.channel+'1.SAC')
            # Horizontal 2 channel
            file2 = datapath / (tstamp+'.'+sta.channel+'2.SAC')
            # Vertical channel
            fileZ = datapath / (tstamp+'.'+sta.channel+'Z.SAC')
            # Pressure channel
            fileP = datapath / (tstamp+'.'+sta.channel+'H.SAC')

            if "P" not in args.channels:

                # If data files exist, continue
                if fileZ.exists() and file1.exists() and file2.exists():
                    if not args.ovr:
                        print(
                            "*   "+tstamp +
                            "*SAC                                 ")
                        print(
                            "*   -> Files already exist, " +
                            "continuing            ")
                        t1 += dt
                        t2 += dt
                        continue

                channels = sta.channel.upper()+'1,'+sta.channel.upper() + \
                    '2,'+sta.channel.upper()+'Z'

                # Get waveforms from client
                try:
                    print("*   "+tstamp +
                          "*SAC                                 ")
                    print("*   -> Downloading Seismic data... ")
                    sth = client.get_waveforms(
                        network=sta.network, station=sta.station,
                        location=sta.location[0], channel=channels,
                        starttime=t1, endtime=t2, attach_response=True)
                    print("*      ...done")
                except:
                    print(" Error: Unable to download ?H? components - " +
                          "continuing")
                    t1 += dt
                    t2 += dt
                    continue

                # Make sure length is ok
                llZ = len(sth.select(component='Z')[0].data)
                ll1 = len(sth.select(component='1')[0].data)
                ll2 = len(sth.select(component='2')[0].data)

                if (llZ != ll1) or (llZ != ll2):
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

            elif "H" not in args.channels:

                # If data files exist, continue
                if fileZ.exists() and fileP.exists():
                    if not args.ovr:
                        print("*   "+tstamp +
                              "*SAC                                 ")
                        print("*   -> Files already exist, " +
                              "continuing            ")
                        t1 += dt
                        t2 += dt
                        continue

                channels = sta.channel.upper() + 'Z'

                # Get waveforms from client
                try:
                    print("*   "+tstamp +
                          "*SAC                                 ")
                    print("*   -> Downloading Seismic data... ")
                    sth = client.get_waveforms(
                        network=sta.network, station=sta.station,
                        location=sta.location[0], channel=channels,
                        starttime=t1, endtime=t2, attach_response=True)
                    print("*      ...done")
                except:
                    print(" Error: Unable to download ?H? components - " +
                          "continuing")
                    t1 += dt
                    t2 += dt
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
                    t1 += dt
                    t2 += dt
                    continue

                # Make sure length is ok
                llZ = len(sth.select(component='Z')[0].data)
                llP = len(stp[0].data)

                if (llZ != llP):
                    print(" Error: lengths not all the same - continuing")
                    t1 += dt
                    t2 += dt
                    continue

                ll = int(dt*stp[0].stats.sampling_rate)

                if np.abs(llZ - ll) > 1:
                    print(" Error: Time series too short - continuing")
                    print(np.abs(llZ - ll))
                    t1 += dt
                    t2 += dt
                    continue

            else:

                # If data files exist, continue
                if (fileZ.exists() and file1.exists() and
                        file2.exists() and fileP.exists()):
                    if not args.ovr:
                        print("*   "+tstamp +
                              "*SAC                                 ")
                        print("*   -> Files already exist, " +
                              "continuing            ")
                        t1 += dt
                        t2 += dt
                        continue

                channels = sta.channel.upper()+'1,'+sta.channel.upper() + \
                    '2,'+sta.channel.upper()+'Z'

                # Get waveforms from client
                try:
                    print("*   "+tstamp +
                          "*SAC                                 ")
                    print("*   -> Downloading Seismic data... ")
                    sth = client.get_waveforms(
                        network=sta.network, station=sta.station,
                        location=sta.location[0], channel=channels,
                        starttime=t1, endtime=t2, attach_response=True)
                    print("*      ...done")
                except:
                    print(" Error: Unable to download ?H? components - " +
                          "continuing")
                    t1 += dt
                    t2 += dt
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
            print("*   -> Removing responses - Seismic data")
            sth.remove_response(pre_filt=args.pre_filt, output='DISP')
            if "P" in args.channels:
                print("*   -> Removing responses - Pressure data")
                stp.remove_response(pre_filt=args.pre_filt)

            # Detrend, filter - seismic data
            sth.detrend('demean')
            sth.detrend('linear')
            sth.filter('lowpass', freq=0.5*args.new_sampling_rate,
                       corners=2, zerophase=True)
            sth.resample(args.new_sampling_rate)

            if "P" in args.channels:
                # Detrend, filter - pressure data
                stp.detrend('demean')
                stp.detrend('linear')
                stp.filter('lowpass', freq=0.5*args.new_sampling_rate,
                           corners=2, zerophase=True)
                stp.resample(args.new_sampling_rate)

            # Extract traces - Z
            trZ = sth.select(component='Z')[0]
            trZ = utils.update_stats(
                trZ, sta.latitude, sta.longitude, sta.elevation, 'Z')
            trZ.write(fileZ, format='sac')

            # Extract traces - H
            if "H" in args.channels:
                tr1 = sth.select(component='1')[0]
                tr2 = sth.select(component='2')[0]
                tr1 = utils.update_stats(
                    tr1, sta.latitude, sta.longitude, sta.elevation, '1')
                tr2 = utils.update_stats(
                    tr2, sta.latitude, sta.longitude, sta.elevation, '2')
                tr1.write(file1, format='sac')
                tr2.write(file2, format='sac')

            # Extract traces - P
            if "P" in args.channels:
                trP = stp[0]
                trP = utils.update_stats(
                    trP, sta.latitude, sta.longitude, sta.elevation, 'P')
                trP.write(fileP, format='sac')

            t1 += dt
            t2 += dt


if __name__ == "__main__":

    # Run main program
    main()
