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

"""

# Import modules and functions
import numpy as np
import os.path
import pickle
import glob
import stdb
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

        # Split into 24-hour long segments
        dt = 3600.*24.

        t1 = tstart
        t2 = tstart + dt

        while t2 <= tend:

            # Time stamp
            tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'

            print(" ")
            print("****************************************************")
            print("* Downloading day-long data for key "+stkey+" and day "+year+"."+jday)

            # Define file names (to check if files already exist)
            file1 = datapath + tstamp + '.' + sta.channel + '1.SAC'
            file2 = datapath + tstamp + '.' + sta.channel + '2.SAC'
            fileZ = datapath + tstamp + '.' + sta.channel + 'Z.SAC'
            fileP = datapath + tstamp + '.' + sta.channel + 'H.SAC'

            # If data file exists, continue
            if glob.glob(fileZ) and glob.glob(file1) and glob.glob(file2) and glob.glob(fileP): 
                if not opts.ovr:
                    print("*   "+tstamp+"*SAC                                 |")
                    print("*   -> Files already exist, continuing            |")
                    t1 += dt
                    t2 += dt
                    continue

            channels = sta.channel.upper() + '1,' + sta.channel.upper() + '2,' + sta.channel.upper() + 'Z'

            # Get waveforms from client
            try:
                print("*   "+tstamp+"*SAC                                 |")
                print("*   -> Downloading Seismic data... ")
                sth = client.get_waveforms(network=sta.network, station=sta.station, location=sta.location[0], \
                        channel=channels, starttime=t1, endtime=t2, attach_response=True)
                print("       ...done")
            except:
                print(" Error: Unable to download ?H? components - continuing")
                t1 += dt
                t2 += dt
                continue
            try:
                print("*   -> Downloading Pressure data...")
                stp = client.get_waveforms(network=sta.network, station=sta.station, location=sta.location[0], \
                        channel='??H', starttime=t1, endtime=t2, attach_response=True)
                print("       ...done")
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
            print("*   -> Removing responses - Seismic data")
            sth.remove_response(pre_filt=opts.pre_filt, output='DISP')
            print("*   -> Removing responses - Pressure data")
            stp.remove_response(pre_filt=opts.pre_filt)

            # Detrend, filter - seismic data
            sth.detrend('demean')
            sth.detrend('linear')
            sth.filter('lowpass', freq=0.5*opts.new_sampling_rate, corners=2, zerophase=True)
            sth.resample(opts.new_sampling_rate)

            # Detrend, filter - pressure data
            stp.detrend('demean')
            stp.detrend('linear')
            stp.filter('lowpass', freq=0.5*opts.new_sampling_rate, corners=2, zerophase=True)
            stp.resample(opts.new_sampling_rate)

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
