#!/usr/bin/env python

# Copyright 2019 Pascal Audet
#
# This file is part of RfPy.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# Import modules and functions
import numpy as np
import os.path
import pickle
import glob
import stdb
from obspy import Stream
from obstools.orient import rf_orient, options
from rfpy import binning


def main():

    # Run Input Parser
    (opts, indb) = options.get_orient_options()

    # Load Database
    db = stdb.io.load_db(fname=indb)

    # Construct station key loop
    allkeys = db.keys()

    # Extract key subset
    if len(opts.stkeys) > 0:
        stkeys = []
        for skey in opts.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = db.keys()

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Define path to see if it exists
        datapath = 'DATA/' + stkey
        if not os.path.isdir(datapath):
            print('Path to '+datapath+' doesn`t exist - continuing')
            continue

        # Temporary print locations
        tlocs = sta.location
        if len(tlocs) == 0:
            tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0:
                tlocs[il] = "--"
        sta.location = tlocs

        # Update Display
        print(" ")
        print(" ")
        print("|===============================================|")
        print("|                   {0:>8s}                    |".format(
            sta.station))
        print("|===============================================|")
        print("|  Station: {0:>2s}.{1:5s}                            |".format(
            sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}  |".format(
            sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(
            sta.longitude, sta.latitude))
        print("|-----------------------------------------------|")

        rfRstream = Stream()
        rfTstream = Stream()

        for folder in os.listdir(datapath):

            filename = datapath+"/"+folder+"/RF_Data.pkl"
            if os.path.isfile(filename):
                file = open(filename, "rb")
                rfdata = pickle.load(file)
                if rfdata[0].stats.snr > opts.snr and \
                        rfdata[0].stats.cc > opts.cc:
                    # if np.std(rfdata[1].data) < 0.2 and \
                    #         np.std(rfdata[2].data) < 0.2:

                    rfRstream.append(rfdata[1])
                    rfTstream.append(rfdata[2])
                file.close()

        if len(rfRstream) == 0:
            continue

        # Remove outliers wrt variance
        # Calculate variance over 30. sec
        nt = int(10./rfRstream[0].stats.delta)
        varR = np.array([np.var(tr.data[0:nt]) for tr in rfRstream])

        # Calculate outliers
        medvarR = np.median(varR)
        madvarR = 1.4826*np.median(np.abs(varR-medvarR))
        robustR = np.abs((varR-medvarR)/madvarR)
        outliersR = np.arange(len(rfRstream))[robustR > 2.5]
        for i in outliersR[::-1]:
            rfRstream.remove(rfRstream[i])
            rfTstream.remove(rfTstream[i])

        # Do the same for transverse
        varT = np.array([np.var(tr.data[0:nt]) for tr in rfTstream])
        medvarT = np.median(varT)
        madvarT = 1.4826*np.median(np.abs(varT-medvarT))
        robustT = np.abs((varT-medvarT)/madvarT)
        outliersT = np.arange(len(rfTstream))[robustT > 2.5]
        for i in outliersT[::-1]:
            rfRstream.remove(rfRstream[i])
            rfTstream.remove(rfTstream[i])

        rfRstream.filter('bandpass', freqmin=opts.fmin,
                         freqmax=opts.fmax, corners=2,
                         zerophase=True)
        rfTstream.filter('bandpass', freqmin=opts.fmin,
                         freqmax=opts.fmax, corners=2,
                         zerophase=True)

        rf_tmp = binning.bin(rfRstream, rfTstream,
                             typ='baz', nbin=opts.nbaz+1)

        azcorr, *_ = rf_orient.decompose(
            rf_tmp[0], rf_tmp[1], opts.t1, opts.t2)
        print("Best fit azcorr: "+
            "{0:5.1f}".format(azcorr))

        azcorr, err_azcorr = rf_orient.get_azcorr(
            rf_tmp[0], rf_tmp[1], opts.t1, opts.t2)
        print("Bootstrap azcorr and uncertainty: "+
            "{0:5.1f}, {1:5.1f}".format(azcorr, err_azcorr))


if __name__ == "__main__":

    # Run main program
    main()
