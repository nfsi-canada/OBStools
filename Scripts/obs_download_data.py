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
from obstools import utils

# Main function
def main():

    dbfile = 'M08A.pkl'

    stationdb = pickle.load(open(dbfile,'rb'))

    sta_key = ['7D.M08A']

    # Get one month of data
    tstart = UTCDateTime('2012-03-01')
    tend = UTCDateTime('2012-03-30')

    # Extract station information from dictionary
    sta = stationdb[sta_key[0]] 

    # Define path to see if it exists
    datapath = 'DATA/' + sta_key[0] + '/'
    if not os.path.isdir(datapath): 
        print('Path to '+datapath+' doesn`t exist - creating it')
        os.makedirs(datapath)

    # Establish client
    client = Client()

    # Split into 24-hour long segments
    dt = 3600.*24.
    new_sampling_rate = 5.

    t1 = tstart
    t2 = tstart + dt

    while t2 <= tend:

        print(t1,t2)
        
        # Time stamp
        tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'
 
        # Define file names (to check if files already exist)
        file1 = datapath + tstamp + '.' + sta.channel + '1.SAC'
        file2 = datapath + tstamp + '.' + sta.channel + '2.SAC'
        fileZ = datapath + tstamp + '.' + sta.channel + 'Z.SAC'
        fileP = datapath + tstamp + '.' + sta.channel + 'H.SAC'

        # If data file exists, continue
        if glob.glob(fileZ) and glob.glob(file1) and glob.glob(file2) and glob.glob(fileP): 
            print('files already exist, continuing')
            t1 += dt
            t2 += dt
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
            t1 += dt
            t2 += dt
            continue
        try:
            print('Downloading Pressure data...')
            stp = client.get_waveforms(network=sta.network, station=sta.station, location=sta.location[0], \
                    channel='??H', starttime=t1, endtime=t2, attach_response=True)
            print('...done')
        except:
            print('not able to download ??H component')
            t1 += dt
            t2 += dt
            continue

        # Make sure length is ok
        llZ = len(sth.select(component='Z')[0].data)
        ll1 = len(sth.select(component='1')[0].data)
        ll2 = len(sth.select(component='2')[0].data)
        llP = len(stp[0].data)

        if (llZ != ll1) or (llZ != ll2) or (llZ != llP):
            print('lengths not all the same - continuing')
            t1 += dt
            t2 += dt
            continue

        ll = int(dt*sth[0].stats.sampling_rate)

        if np.abs(llZ - ll) > 1:
            print('Time series too short - continuing')
            print(np.abs(llZ - ll))
            t1 += dt
            t2 += dt
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
