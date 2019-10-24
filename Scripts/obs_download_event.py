#!/usr/bin/env python

'''
PROGRAM request_IRIS.py

Program to request day-long seismograms from the IRIS archive.

'''

# Import modules and functions
import numpy as np
import os.path
import pickle
import glob
from obspy import UTCDateTime
from obspy.core import AttribDict
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth as epi
from obspy.geodetics import kilometer2degrees as k2d
from obstools import utils
from obstools import EventStream

# Main function
def main():

    dbfile = 'M08A.pkl'

    stationdb = pickle.load(open(dbfile,'rb'))

    # Initialize Taup Model
    tpmodel = TauPyModel(model='iasp91')

    sta_key = ['7D.M08A']

    # Get one month of data
    tstart = UTCDateTime('2012-03-09')
    tend = UTCDateTime('2012-03-10')

    # Extract station information from dictionary
    sta = stationdb[sta_key[0]] 

    # Define path to see if it exists
    eventpath = 'EVENTS/' + sta_key[0] + '/'
    if not os.path.isdir(eventpath): 
        print('Path to '+eventpath+' doesn`t exist - creating it')
        os.makedirs(eventpath)

    # Establish client
    client = Client()

    # Get catalogue using deployment start and end
    cat = client.get_events(starttime=tstart, endtime=tend,    \
                minmagnitude=6.0, maxmagnitude=7.0)

    # Total number of events in Catalogue
    nevtT = len(cat)
    print("|  Found {0:5d} possible events                  |".format(nevtT))
    ievs = range(0, nevtT)

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

        # Get Travel times (Careful: here dep is in meters)
        tt = tpmodel.get_travel_times(distance_in_degree=gac, \
            source_depth_in_km=dep/1000.)

        # Loop over all times in tt
        for t in tt:

            # Extract SKS arrival
            if t.name == 'P':

                # Add SKS phase to Split object
                tarrival = t.time

                # Break out of loop 
                break

        t1 = time # + tarrival - 800.
        t2 = t1 + window

        print(t1,t2)
        
        # Time stamp
        tstamp = str(time.year).zfill(4)+'.'+str(time.julday).zfill(3)+'.'
        tstamp = tstamp + str(time.hour).zfill(2)+'.'+str(time.minute).zfill(2)
 
        # Define file names (to check if files already exist)
        filename = eventpath + tstamp + '.event.pkl'

        # If data file exists, continue
        if glob.glob(filename): 
            print('files already exist, continuing')
            continue

        channels = sta.channel.upper() + '1,' + sta.channel.upper() + '2,' + sta.channel.upper() + 'Z'

        # Get waveforms from client
        try:
            print('Downloading Seismic data...')
            sth = client.get_waveforms(network=sta.network, station=sta.station, location=sta.location[0], \
                    channel=channels, starttime=t1, endtime=t2, attach_response=True)
            print('...done')
        except:
            print('not able to download BH? components')
            continue
        try:
            print('Downloading Pressure data...')
            stp = client.get_waveforms(network=sta.network, station=sta.station, location=sta.location[0], \
                    channel='??H', starttime=t1, endtime=t2, attach_response=True)
            print('...done')
        except:
            print('not able to download ??H component')
            continue

        # Make sure length is ok
        llZ = len(sth.select(component='Z')[0].data)
        ll1 = len(sth.select(component='1')[0].data)
        ll2 = len(sth.select(component='2')[0].data)
        llP = len(stp[0].data)

        if (llZ != ll1) or (llZ != ll2) or (llZ != llP):
            print('lengths not all the same - continuing')
            continue

        ll = int(window*sth[0].stats.sampling_rate)

        if np.abs(llZ - ll) > 1:
            print('Time series too short - continuing')
            print(np.abs(llZ - ll))
            continue

        # Remove responses
        print('Removing responses - Seismic data')
        sth.remove_response(pre_filt=[0.001,0.005, 45., 50.], output='DISP')
        print('Removing responses - Pressure data')
        stp.remove_response(pre_filt=[0.001,0.005, 45., 50.])

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

        eventstream = EventStream(sta, sth, stp, tstamp, lat, lon, time, window, new_sampling_rate)

        eventstream.save(filename)


if __name__ == "__main__":

    # Run main program
    main()
