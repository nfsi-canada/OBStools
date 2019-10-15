# Import modules and functions
import os
import fnmatch
import numpy as np
from obspy.core import read, Stream, Trace, AttribDict
from obspy import UTCDateTime
import pickle
from obstools import QC as qc

def main():

    dbfile = 'M08A.pkl'

    stationdb = pickle.load(open(dbfile,'rb'))

    sta_key = ['7D.M08A']

    tstart = UTCDateTime('2012-03-01')
    tend = UTCDateTime('2012-03-10')

    # Extract station information from dictionary
    sta = stationdb[sta_key[0]] 

    datapath = 'DATA/' + sta.network + '/' + sta.station + '/'

    # Get all components
    trN1, trN2, trNZ, trNP = get_data(datapath, tstart, tend)

    # Window size 
    window = 7200.
    overlap = 0.3

    # minimum numer of windows
    minwin = 10

    # NFFT
    stats = trN1[0].stats

    # Time axis
    taxis = np.arange(0., window, trN1[0].stats.delta)

    for tr1, tr2, trZ, trP in zip(trN1, trN2, trNZ, trNP):

        # Quality control to identify outliers
        goodwins = qc.QC_daily_spectra(tr1, tr2, trZ, trP, window, overlap, fig_QC_specs=False)

        # Check if we have enough good windows
        nwin = np.sum(goodwins)
        if nwin < minwin:
            print('Too few good data segments to calculate day spectra')
            # continue
        else:
            print('{0} good windows. Proceeding...'.format(nwin))

        specprop = qc.spec_calc_daily_spectra( \
            tr1, tr2, trZ, trP, goodwins, window, overlap, fig_av_specs=False)



# Function to read all available receiver functions that meet SNR threshold
def get_data(datapath, tstart, tend):

    t1 = tstart
    t2 = tend

    # Define empty streams
    trN1 = Stream()
    trN2 = Stream()
    trNZ = Stream()
    trNP = Stream()

    while t1 < t2:

        tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'

        # Loop through directory and load files
        for file in os.listdir(datapath):
            if fnmatch.fnmatch(file, '*' + tstamp + '*1.SAC'):
                tr = read(datapath + file)
                trN1.append(tr[0])
            elif fnmatch.fnmatch(file, '*' + tstamp + '*2.SAC'):
                tr = read(datapath + file)
                trN2.append(tr[0])
            elif fnmatch.fnmatch(file, '*' + tstamp + '*Z.SAC'):
                tr = read(datapath + file)
                trNZ.append(tr[0])
            elif fnmatch.fnmatch(file, '*' + tstamp + '*H.SAC'):
                tr = read(datapath + file)
                trNP.append(tr[0])

        t1 += 3600.*24.

    return trN1, trN2, trNZ, trNP


if __name__ == "__main__":

    # Run main program
    main()
