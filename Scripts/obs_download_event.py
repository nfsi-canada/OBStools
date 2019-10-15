'''
PROGRAM obs_get_data.py

Grabs three component data from the IRIS archive for OBS stations

Station selection is specified by a network and station codes.
The data base is provided in stations_db.py as a dictionary.

'''

# Import modules and functions
import numpy as np
import numpy.fft as ft
import matplotlib.pyplot as plt
import obspy.core
import os.path
import pickle
from math import pi, sin
from obspy.core.event import readEvents
from obspy import UTCDateTime
from obspy.core import read, Stream, Trace
from obspy.fdsn import Client
from obspy.core.util.geodetics import gps2DistAzimuth as epi
from obspy.core.util.geodetics import kilometer2degrees as k2d
from obspy.taup.taup import getTravelTimes as gtt
from obspy.signal.rotate import rotate_NE_RT 
from obs import obs_proc as obsp

# Plot intermediate steps
plot = False

# Main function
def process(db, sta_key):

    # Extract station information from dictionary
    sta = db[sta_key]

    # Define path to see if it exists
    trfpath = 'OBS_DATA/'+sta.station
    if not os.path.isdir(trfpath): 
        print 'Path to '+trfpath+' doesn`t exist - creating it'
        os.makedirs(trfpath)

    # Establish client
    client = Client()

    # Get catalogue using deployment start and end
    cat = client.get_events(starttime=sta.dstart, endtime=sta.dend, minmagnitude=6.0)

    # Read catalogue
    for ev in cat:

        # Extract time, coordinates and depth of events
        time = ev.origins[0].time
        lat = ev.origins[0].latitude
        lon = ev.origins[0].longitude
        dep = ev.origins[0].depth
        
        # Define time stamp
        yr = str(time.year).zfill(4)
        jd = str(time.julday).zfill(3)
        hr = str(time.hour).zfill(2)
        tstamp = yr+jd+hr

        # Define file names (to check if files already exist)
        fileN1 = trfpath+'/OBS_N1.'+sta.network+'.'+sta.station+'.'+tstamp+'.'+'sac'
        fileN2 = trfpath+'/OBS_N2.'+sta.network+'.'+sta.station+'.'+tstamp+'.'+'sac'
        fileNZ = trfpath+'/OBS_NZ.'+sta.network+'.'+sta.station+'.'+tstamp+'.'+'sac'
        fileNP = trfpath+'/OBS_NP.'+sta.network+'.'+sta.station+'.'+tstamp+'.'+'sac'
        file1 = trfpath+'/OBS1.'+sta.network+'.'+sta.station+'.'+tstamp+'.'+'sac'
        file2 = trfpath+'/OBS2.'+sta.network+'.'+sta.station+'.'+tstamp+'.'+'sac'
        fileZ = trfpath+'/OBSZ.'+sta.network+'.'+sta.station+'.'+tstamp+'.'+'sac'
        fileP = trfpath+'/OBSP.'+sta.network+'.'+sta.station+'.'+tstamp+'.'+'sac'
        
        # If RF file exists, continue
        if os.path.isfile(fileZ): continue

        # Calculate epicentral distance
        epi_dist, az, baz = epi(lat, lon, sta.stla, sta.stlo)
        epi_dist /= 1000
        gac = k2d(epi_dist)
 
        # If distance between 30 and 90 deg:
        if gac>30. and gac<90.:

            #print
            #print 'EARTHQUAKE INFO'
            print
            print 'Station:     ', sta.station
            print 'Time:        ', time
            print 'Backazimuth: ', baz
            print 'Event depth: ', dep/1000.
            if not dep:
                dep = 10000.

            # Get Travel times (Careful: here dep is in meters)
            tt = gtt(delta=gac, depth=dep/1000., model='iasp91')

            # Loop over all times in tt
            for t in tt:

                # Extract time of P arrival
                if t['phase_name']=='P':
                    
                    # Extract time, phase name and slowness
                    tp = t['time']
                    ph = t['phase_name']
                    slow = t['dT/dD']/111.

                    # Break out of loop 
                    break

            # Define start and end times for requests
            tstart = time+tp-6.*3600.
            tend = time+tp
            
            # Get waveforms from client
            try:
                sth = client.get_waveforms(sta.network,sta.station,'*', \
                        sta.cha+'?',tstart, tend, attach_response=True)
            except Exception as e:
                continue
            try:
                stp = client.get_waveforms(sta.network,sta.station,'*', \
                        '?XH',tstart, tend, attach_response=True)
            except Exception as e:
                continue
            sth.remove_response(output='DISP')
            stp.remove_response()

            # Detrend, filter
            sth.detrend('demean')
            sth.detrend('linear')
            sth.filter('lowpass', freq=0.5*2.0, corners=2, zerophase=True)
            sth.resample(2.0)

            stp.detrend('demean')
            stp.detrend('linear')
            stp.filter('lowpass', freq=0.5*2.0, corners=2, zerophase=True)
            stp.resample(2.0)

            # Extract traces
            trZ = sth.select(component='Z')[0]
            tr1 = sth.select(component='1')[0]
            tr2 = sth.select(component='2')[0]
            trP = stp.select(channel='?XH')[0]

            # Window size 
            ws = int(2000./tr1.stats.delta)

            # Step size
            ss = int(2000./tr1.stats.delta)

            # Update stats
            tr1 = update_stats(tr1, sta.stla, sta.stlo, sta.stel, ws, ss)
            tr2 = update_stats(tr2, sta.stla, sta.stlo, sta.stel, ws, ss)
            trZ = update_stats(trZ, sta.stla, sta.stlo, sta.stel, ws, ss)
            trP = update_stats(trP, sta.stla, sta.stlo, sta.stel, ws, ss)

            # Save traces
            tr1.write(fileN1,format='sac')
            tr2.write(fileN2,format='sac')
            trZ.write(fileNZ,format='sac')
            trP.write(fileNP,format='sac')

            # Define start and end times for requests
            tstart = time+tp-120.
            tend = time+tp+120.
            
            # Get waveforms from client
            try:
                sth = client.get_waveforms(sta.network,sta.station,'*', \
                        sta.cha+'?',tstart, tend, attach_response=True)
            except Exception as e:
                continue
            try:
                stp = client.get_waveforms(sta.network,sta.station,'*', \
                        '?XH',tstart, tend, attach_response=True)
            except Exception as e:
                continue

            sth.remove_response(output='DISP')
            stp.remove_response()

            # Detrend, filter
            sth.detrend('demean')
            sth.detrend('linear')
            sth.filter('lowpass', freq=0.5*2.0, corners=2, zerophase=True)
            sth.resample(2.0)

            stp.detrend('demean')
            stp.detrend('linear')
            stp.filter('lowpass', freq=0.5*2.0, corners=2, zerophase=True)
            stp.resample(2.0)

            # Extract traces
            trZ = sth.select(component='Z')[0]
            tr1 = sth.select(component='1')[0]
            tr2 = sth.select(component='2')[0]
            trP = stp.select(channel='?XH')[0]

            # Update stats
            tr1 = update_stats(tr1, sta.stla, sta.stlo, sta.stel, ws, ss)
            tr2 = update_stats(tr2, sta.stla, sta.stlo, sta.stel, ws, ss)
            trZ = update_stats(trZ, sta.stla, sta.stlo, sta.stel, ws, ss)
            trP = update_stats(trP, sta.stla, sta.stlo, sta.stel, ws, ss)

            # Save traces
            tr1.write(file1,format='sac')
            tr2.write(file2,format='sac')
            trZ.write(fileZ,format='sac')
            trP.write(fileP,format='sac')


def rotate_12_NE(tr1, tr2, az):

    a = az*np.pi/180.
    rot_mat = np.array([[np.cos(a), np.sin(a)],
        [-np.sin(a), np.cos(a)]])

    v12 = np.array([tr2,tr1])
    vne = np.dot(rot_mat, v12)
    trE = vne[0,:]
    trN = vne[1,:]
    
    return trN, trE


def update_stats(tr, stla, stlo, stel, ws, ss):
    tr.stats.sac = obspy.core.AttribDict()
    tr.stats.sac.stla = stla
    tr.stats.sac.stlo = stlo
    tr.stats.sac.stel = stel
    tr.stats.sac.user8 = ws
    tr.stats.sac.user9 = ss
    return tr


# Loads station db and builds attribute dict of station stats
def load_db(fname):
    db = pickle.load(open(fname, 'rb'))
    for k, v in db.items():
        db[k] = meta_data(v)
    return db


# Attribute dict class
class meta_data(dict):
    def __init__(self, stats):
        self.__dict__ = self
        self.network = stats[0]
        self.station = stats[1]
        self.stla = stats[2]
        self.stlo = stats[3]
        self.stel = stats[4]
        self.azim = stats[5]
        self.cha = stats[6]
        self.dstart = stats[7]
        self.dend = stats[8]


###############################
# Choose one station to process

dbfile = 'stations_obs.pkl'
stationdb = load_db(dbfile)

sta_keys = ['J39C']

for key in sta_keys:
    process(stationdb, key)

###################
