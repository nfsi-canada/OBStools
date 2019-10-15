'''
PROGRAM obs_prep.py

'''

# Import modules and functions
import os
import fnmatch
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import read, Stream, Trace, AttribDict
import pickle
from obs import obs_proc as obsp
from obs import obs_utils as obsu


def get_obs(db, sta_key, tstart, tend):

    # Extract station information from dictionary
    sta = db[sta_key] 

    # Get all components
    trN1, trN2, trNZ, trNP, tr1, tr2, trZ, trP = \
            get_data('OBS_DATA/'+sta.station+'/',sta.network,sta.station, \
                tstart, tend)
   
    # Prepare data
    for i in range(len(trN1)):

        # Get tilt direction
        obsu.tilt_direction(trN1[i], trN2[i], trNZ[i], trNP[i])

        # # Transfer function corrections for tilt
        # trNZ_p, trZ_p = obsp.obs_proc_tilt(trN1[i], trN2[i], trNZ[i],\
        #         tr1[i], tr2[i], trZ[i])

        # # Transfer function corrections for compliance
        # trZ_f = obsp.obs_proc_compliance(trNZ_p, trNP[i], trZ_p, trP[i])
        # #trZ_f = obsp.obs_proc_compliance(trNZ[i], trNP[i], trZ[i], trP[i])

    return 


# Function to read all available receiver functions that meet SNR threshold
def get_data(filedir, net, sta, tstart, tend):

    # Define empty streams
    trN1 = Stream()
    trN2 = Stream()
    trNZ = Stream()
    trNP = Stream()
    tr1 = Stream()
    tr2 = Stream()
    trZ = Stream()
    trP = Stream()

    # Loop through directory and load files
    for file in os.listdir(filedir):
        if fnmatch.fnmatch(file, 'OBS_N1.'+net+'.'+sta+'.*.sac'):
            tr = read(filedir+file)
            trN1.append(tr[0])
        elif fnmatch.fnmatch(file, 'OBS_N2.'+net+'.'+sta+'.*.sac'):
            tr = read(filedir+file)
            trN2.append(tr[0])
        elif fnmatch.fnmatch(file, 'OBS_NZ.'+net+'.'+sta+'.*.sac'):
            tr = read(filedir+file)
            trNZ.append(tr[0])
        elif fnmatch.fnmatch(file, 'OBS_NP.'+net+'.'+sta+'.*.sac'):
            tr = read(filedir+file)
            trNP.append(tr[0])
        elif fnmatch.fnmatch(file, 'OBS1.'+net+'.'+sta+'.*.sac'):
            tr = read(filedir+file)
            tr1.append(tr[0])
        elif fnmatch.fnmatch(file, 'OBS2.'+net+'.'+sta+'.*.sac'):
            tr = read(filedir+file)
            tr2.append(tr[0])
        elif fnmatch.fnmatch(file, 'OBSZ.'+net+'.'+sta+'.*.sac'):
            tr = read(filedir+file)
            trZ.append(tr[0])
        elif fnmatch.fnmatch(file, 'OBSP.'+net+'.'+sta+'.*.sac'):
            tr = read(filedir+file)
            trP.append(tr[0])

    return trN1, trN2, trNZ, trNP, tr1, tr2, trZ, trP


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


#######################
# Execute main program

dbfile = 'stations.pkl'
stationdb = load_db(dbfile)

sta_keys = ['G03A']
tstart = UTCDateTime('2011-12-01')
tend = UTCDateTime('2011-12-01')

for key in sta_keys:

    print('--------------------------------------------')
    print('Pre-processing OBS data '+key)

    # Call function to calculate H-k solution
    get_obs(stationdb, key, tstart, tend)


