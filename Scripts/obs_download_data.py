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

# Main function
def process(db, sta_key):

    sta = db[sta_key]

    # Define path to see if it exists
    datapath = 'DATA/' + sta.network + '/' + sta.station + '/'
    if not os.path.isdir(datapath): 
        print('Path to '+datapath+' doesn`t exist - creating it')
        os.makedirs(datapath)

    # Establish client
    client = Client()

    # Get one month of data
    tstart = UTCDateTime('2012-03-01')
    tend = UTCDateTime('2012-03-30')

    # Split into 24-hour long segments
    dt = 3600.*24.
    new_sampling_rate = 5.

    t1 = tstart
    t2 = tstart + dt

    print(sta.station, sta.network, sta.channel, sta.location)

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
        if glob.glob(fileZ) and glob.glog(file1) and glob.glob(file2) and glob.glob(fileP): 
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
        llZ = len(sth.select(component='Z')[0].data) < int(dt*sth[0].stats.sampling_rate)
        ll1 = len(sth.select(component='1')[0].data) < int(dt*sth[0].stats.sampling_rate)
        ll2 = len(sth.select(component='2')[0].data) < int(dt*sth[0].stats.sampling_rate)
        llP = len(stp[0].data) < int(dt*stp[0].stats.sampling_rate)

        if llZ or ll1 or ll2 or llP:
            print('Time series too short - continuing')
            t1 += dt
            t2 += dt
            continue

        # Remove responses
        print('Removing responses - Seismic data')
        sth.remove_response(output='DISP')
        print('Removing responses - Pressure data')
        stp.remove_response()

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
        tr1 = update_stats(tr1, sta.latitude, sta.longitude, sta.elevation)
        tr2 = update_stats(tr2, sta.latitude, sta.longitude, sta.elevation)
        trZ = update_stats(trZ, sta.latitude, sta.longitude, sta.elevation)
        trP = update_stats(trP, sta.latitude, sta.longitude, sta.elevation)

        # Save traces
        tr1.write(file1, format='sac')
        tr2.write(file2, format='sac')
        trZ.write(fileZ, format='sac')
        trP.write(fileP, format='sac')

        t1 += dt
        t2 += dt

def update_stats(tr, stla, stlo, stel):
    tr.stats.sac = AttribDict()
    tr.stats.sac.stla = stla
    tr.stats.sac.stlo = stlo
    tr.stats.sac.stel = stel
    return tr


###############################
# Choose one station to process

dbfile = 'M08A.pkl'

stationdb = pickle.load(open(dbfile,'rb'))

#sta_keys = ['WHY', 'INK', 'EPYK']
#sta_keys = ['F55A', 'G55A', 'G57A', 'H56A']
sta_keys = ['7D.M08A']

for key in sta_keys:
    process(stationdb, key)

###################
