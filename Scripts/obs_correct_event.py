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

    # Path where transfer functions are located
    transpath = 'TF_STA/' + sta_key[0] + '/'

    # Path where event data are located
    eventpath = 'EVENTS/' + sta_key[0] + '/'

    # Find all files in directories
    event_files = os.listdir(eventpath)
    trans_files = os.listdir(transpath)

    # # List of possible transfer functions for Daily files
    # TF_list_day = {'ZP': True, 'Z1':False, 'Z2-1':False, 'ZP-21':False, 'ZH':True, 'ZP-H':False}

    # Cycle through available files
    for eventfile in event_files:

        if eventfile[0]=='.':
            continue

        evprefix = eventfile.split('.')
        evstamp = evprefix[0]+'.'+evprefix[1]+'.'

        # Load event file
        file = open(eventpath+eventfile, 'rb')
        eventstream = pickle.load(file)
        file.close()

        # if fig_event_raw:
        plot.fig_event_raw(eventstream)

        # Cycle through corresponding TF files
        for transfile in trans_files:

            if transfile[0]=='.':
                continue

            tfprefix = transfile.split('transfunc')[0]
            if len(tfprefix) > 9:
                yr1 = tfprefix.split('-')[0].split('.')[0]
                jd1 = tfprefix.split('-')[0].split('.')[1]
                yr2 = tfprefix.split('-')[1].split('.')[0]
                jd2 = tfprefix.split('-')[1].split('.')[1]
                if evprefix[0]>=yr1 and evprefix[0] <=yr2:
                    if evprefix[1]>=jd1 and evprefix[1] <=jd2:
                        file = open(transpath+transfile, 'rb')
                        print(transpath+transfile)
                        tfaverage = pickle.load(file)
                        file.close()

                        # List of possible transfer functions for station average files
                        TF_list = {'ZP': True, 'Z1':True, 'Z2-1':True, 'ZP-21':True, 'ZH':False, 'ZP-H':False}
                        eventstream.correct_data(tfaverage, TF_list)

                        correct = eventstream.correct
                        plot.fig_event_corrected(eventstream, TF_list)

            else:
                if tfprefix==evstamp:
                    file = open(transpath+transfile, 'rb')
                    print(transpath+transfile)
                    tfaverage = pickle.load(file)
                    file.close()

                    # List of possible transfer functions for station average files
                    TF_list = {'ZP': True, 'Z1':True, 'Z2-1':True, 'ZP-21':True, 'ZH':True, 'ZP-H':True}
                    eventstream.correct_data(tfaverage, TF_list)

                    correct = eventstream.correct
                    plot.fig_event_corrected(eventstream, TF_list)


if __name__ == "__main__":

    # Run main program
    main()

