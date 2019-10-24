#!/usr/bin/env python

# Import modules and functions
import os
import numpy as np
from obspy import UTCDateTime
import pickle
from obstools import StaNoise, Power, Cross, Rotation, TFNoise
from obstools import utils, plot

def main():

    dbfile = 'M08A.pkl'

    stationdb = pickle.load(open(dbfile,'rb'))

    sta_key = ['7D.M08A']

    tstart = UTCDateTime('2012-03-01')
    tend = UTCDateTime('2012-03-30')

    # Extract station information from dictionary
    sta = stationdb[sta_key[0]] 

    # Path where spectra are located
    specpath = 'SPECTRA/' + sta_key[0] + '/'

    # Path where average spectra are located
    avstpath = 'AVG_STA/' + sta_key[0] + '/'

    # Path where transfer functions will be located
    tfpath = 'TF_STA/' + sta_key[0] + '/'
    if not os.path.isdir(tfpath): 
        print('Path to '+tfpath+' doesn`t exist - creating it')
        os.makedirs(tfpath)

    # Filename for output transfer functions
    dstart = str(tstart.year).zfill(4)+'.'+str(tstart.julday).zfill(3)+'-'
    dend = str(tend.year).zfill(4)+'.'+str(tend.julday).zfill(3)+'.'
    fileavst = avstpath + dstart + dend + 'avg_sta.pkl'


    # Find all files in directories
    spectra_files = os.listdir(specpath)
    average_files = os.listdir(avstpath)

    # List of possible transfer functions for Daily files
    TF_list_day = {'ZP': True, 'Z1':True, 'Z2-1':True, 'ZP-21':True, 'ZH':True, 'ZP-H':True}

    day_transfer_functions = []

    # Cycle through available files
    for filespec in spectra_files:

        year = filespec.split('.')[0]
        jday = filespec.split('.')[1]

        print('Calculating transfer functions for key '+sta_key[0]+' and day '+year+'.'+jday)
        tstamp = year+'.'+jday+'.'
        filename = tfpath + tstamp + 'transfunc.pkl'

        # Load file
        file = open(specpath+filespec, 'rb')
        daynoise = pickle.load(file)
        file.close()

        # Load spectra into TFNoise object
        daytransfer = TFNoise(daynoise.f, daynoise.power, daynoise.cross, daynoise.rotation, TF_list_day)

        # Calculate the transfer functions
        daytransfer.transfer_func()

        # Store the frequency axis
        f = daytransfer.f

        # Append to list of transfer functions
        day_transfer_functions.append(daytransfer.transfunc)

        # Save daily transfer functions to file
        daytransfer.save(filename)

    # # Convert to numpy ndarray
    # day_transfer_functions = np.array(day_transfer_functions)

    # print(day_transfer_functions[0]['ZP']['TF_ZP'])
    # print(day_transfer_functions)

    # List of possible transfer functions for station average files
    TF_list_sta = {'ZP': True, 'Z1':True, 'Z2-1':True, 'ZP-21':True, 'ZH':False, 'ZP-H':False}

    # Cycle through available files
    for fileavst in average_files:

        name = fileavst.split('avg_sta')

        print('Calculating transfer functions for key '+sta_key[0]+' and range '+name[0])
        filename = tfpath + name[0] + 'transfunc.pkl'

        # Load file
        file = open(avstpath+fileavst, 'rb')
        stanoise = pickle.load(file)
        file.close()

        # Load spectra into TFNoise object - no Rotation object for station averages
        rotation = Rotation(None, None, None)
        statransfer = TFNoise(stanoise.f, stanoise.power, stanoise.cross, rotation, TF_list_sta)

        # Calculate the transfer functions
        statransfer.transfer_func()

        # Store the frequency axis
        f = statransfer.f

        # Extract the transfer functions
        sta_transfer_functions = statransfer.transfunc

        # Save average transfer functions to file
        statransfer.save(filename)

    # if fig_TF:
    plot.fig_TF(f, day_transfer_functions, TF_list_day, sta_transfer_functions, TF_list_sta, key=sta_key[0])


if __name__ == "__main__":

    # Run main program
    main()

