#!/usr/bin/env python

# Import modules and functions
import os
import numpy as np
from obspy import UTCDateTime
import pickle
from obstools import DayNoise
from obstools import utils

def main():

    dbfile = 'M08A.pkl'

    stationdb = pickle.load(open(dbfile,'rb'))

    sta_key = ['7D.M08A']

    tstart = UTCDateTime('2012-03-01')
    tend = UTCDateTime('2012-03-31')

    # Extract station information from dictionary
    sta = stationdb[sta_key[0]] 

    # Path where data are located
    datapath = 'DATA/' + sta_key[0] + '/'

    # Path where spectra will be saved
    specpath = 'SPECTRA/' + sta_key[0] + '/'
    if not os.path.isdir(specpath): 
        print('Path to '+specpath+' doesn`t exist - creating it')
        os.makedirs(specpath)

    # Get all components
    trN1, trN2, trNZ, trNP = utils.get_data(datapath, tstart, tend)

    # Window size 
    window = 7200.
    overlap = 0.3

    # minimum numer of windows
    minwin = 10

    # NFFT
    stats = trN1[0].stats

    # Time axis
    taxis = np.arange(0., window, trN1[0].stats.delta)

    # Cycle through available data
    for tr1, tr2, trZ, trP in zip(trN1, trN2, trNZ, trNP):

        year = str(tr1.stats.starttime.year).zfill(4)
        jday = str(tr1.stats.starttime.julday).zfill(3)

        print('Calculating noise spectra for key '+sta_key[0]+' and day '+year+'.'+jday)
        tstamp = year+'.'+jday+'.'
        filename = specpath + tstamp + 'spectra.pkl'

        if os.path.exists(filename):
            print('file '+filename+' exists - continuing')
            # continue

        # Initialize instance of DayNoise
        daynoise = DayNoise(tr1, tr2, trZ, trP, window, overlap, key=sta_key[0])

        # Quality control to identify outliers
        daynoise.QC_daily_spectra(fig_QC=False, debug=False)

        # Check if we have enough good windows
        nwin = np.sum(daynoise.goodwins)
        if nwin < minwin:
            print('Too few good data segments to calculate day spectra')
            # continue
        else:
            print('{0} good windows. Proceeding...'.format(nwin))

        # Average spectra for good windows
        daynoise.average_daily_spectra(fig_average=False, fig_coh_ph=False)

        # Save to file
        daynoise.save(filename)



if __name__ == "__main__":

    # Run main program
    main()
