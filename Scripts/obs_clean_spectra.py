#!/usr/bin/env python

# Copyright 2019 Pascal Audet & Helen Janiszewski
#
# This file is part of OBStools.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Program obs_clean_spectra.py
----------------------------

Extracts daily spectra calculated from `obs_daily_spectra.py` and 
flags days for which the daily averages are outliers from the PSD properties. 
Further averages the spectra over the whole period specified by `--start`
and `--end`.

Station selection is specified by a network and 
station code. The data base is provided as a 
`StDb` dictionary.

Usage
-----

.. code-block::

    $ obs_clean_spectra.py -h
    Usage: obs_clean_spectra.py [options] <station database>

    Script used to extract daily spectra calculated from `obs_daily_spectra.py`
    and flag days for outlier PSDs and calculate spectral averages of the
    corresponding Fourier transforms over the entire time period specified. The
    stations are processed one by one and the data are stored to disk.

    Options:
      -h, --help         show this help message and exit
      --keys=STKEYS      Specify a comma separated list of station keys for which
                         to perform the analysis. These must be contained within
                         the station database. Partial keys will be used to match
                         against those in the dictionary. For instance, providing
                         IU will match with all stations in the IU network.
                         [Default processes all stations in the database]
      -O, --overwrite    Force the overwriting of pre-existing data. [Default
                         False]

      Parameter Settings:
        Miscellaneous default values and settings

        --freq-band=PD   Specify comma-separated frequency limits (float, in Hz)
                         over which to calculate spectral features used in
                         flagging the days/windows. [Default 0.004,2.0]
        --tolerance=TOL  Specify parameter for tolerance threshold. If spectrum >
                         std*tol, window is flagged as bad. [Default 1.5]
        --alpha=ALPHA    Confidence level for f-test, for iterative flagging of
                         windows. [Default 0.05, or 95% confidence]

      Figure Settings:
        Flags for plotting figures

        --figQC          Plot Quality-Control figure. [Default does not plot
                         figure]
        --debug          Plot intermediate steps for debugging. [Default does not
                         plot figure]
        --figAverage     Plot daily average figure. [Default does not plot figure]
        --figCoh         Plot Coherence and Phase figure. [Default does not plot
                         figure]
        --figCross       Plot cross-spectra figure. [Default does not plot figure]

      Time Search Settings:
        Time settings associated with searching for day-long seismograms

        --start=STARTT   Specify a UTCDateTime compatible string representing the
                         start day for the data search. This will override any
                         station start times. [Default start date of each station
                         in database]
        --end=ENDT       Specify a UTCDateTime compatible string representing the
                         start time for the event search. This will override any
                         station end times. [Default end date of each station in
                         database]

"""

# Import modules and functions
import os
import sys
import numpy as np
import pickle
import stdb
from obstools import StaNoise, Power, Cross, Rotation
from obstools import utils, options, plot

def main():

    # Run Input Parser
    (opts, indb) = options.get_cleanspec_options()

    # Load Database
    db = stdb.io.load_db(fname=indb)

    # Construct station key loop
    allkeys = db.keys()
    sorted(allkeys)

    # Extract key subset
    if len(opts.stkeys) > 0:
        stkeys = []
        for skey in opts.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = db.keys()
        sorted(stkeys)

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Path where spectra are located
        specpath = 'SPECTRA/' + stkey + '/'
        if not os.path.isdir(avstpath): 
            print("Path to '+specpath+' doesn`t exist - aborting")
            sys.exit()

        # Path where average spectra will be saved
        avstpath = 'AVG_STA/' + stkey + '/'
        if not os.path.isdir(avstpath): 
            print("Path to '+avstpath+' doesn`t exist - creating it")
            os.makedirs(avstpath)

        # Get catalogue search start time
        if opts.startT is None:
            tstart = sta.startdate
        else:
            tstart = opts.startT

        # Get catalogue search end time
        if opts.endT is None:
            tend = sta.startdate
        else:
            tend = opts.endT

        if tstart > sta.enddate or tend < sta.startdate:
            continue

        # Temporary print locations
        tlocs = sta.location
        if len(tlocs) == 0: tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0: tlocs[il] = "--"
        sta.location = tlocs

        # Update Display
        print(" ")
        print(" ")
        print("|===============================================|")
        print("|===============================================|")
        print("|                   {0:>8s}                    |".format(sta.station))
        print("|===============================================|")
        print("|===============================================|")
        print("|  Station: {0:>2s}.{1:5s}                            |".format(sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}  |".format(sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(sta.longitude, sta.latitude))
        print("|      Start time: {0:19s}          |".format(sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|      End time:   {0:19s}          |".format(sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|-----------------------------------------------|")


        # Filename for output average spectra
        dstart = str(tstart.year).zfill(4)+'.'+str(tstart.julday).zfill(3)+'-'
        dend = str(tend.year).zfill(4)+'.'+str(tend.julday).zfill(3)+'.'
        fileavst = avstpath + dstart + dend + 'avg_sta.pkl'

        # Containers for power and cross spectra
        c11_all = []; c22_all = []; cZZ_all = []; cPP_all = []
        c12_all = []; c1Z_all = []; c1P_all = []; c2Z_all = []; c2P_all = []; cZP_all = []
        cHH_all = []; cHZ_all = []; cHP_all = []
        coh_all = []; ph_all = []
        coh_12_all = []; coh_1Z_all = []; coh_1P_all = []; coh_2Z_all = []; coh_2P_all = []; coh_ZP_all = []
        ph_12_all = []; ph_1Z_all = []; ph_1P_all = []; ph_2Z_all = []; ph_2P_all = []; ph_ZP_all = []
        ad_12_all = []; ad_1Z_all = []; ad_1P_all = []; ad_2Z_all = []; ad_2P_all = []; ad_ZP_all = []
        nwins = []

        t1 = tstart

        # Loop through each day withing time range
        while t1 < tend:

            year = str(t1.year).zfill(4)
            jday = str(t1.julday).zfill(3)

            print('Calculating noise spectra for key '+stkey+' and day '+year+'.'+jday)
            tstamp = year+'.'+jday+'.'
            filespec = specpath + tstamp + 'spectra.pkl'

            # Load file if it exists
            if os.path.exists(filespec):
                print('file '+filespec+' exists - loading')
                file = open(filespec, 'rb')
                daynoise = pickle.load(file)
                file.close()
            else:
                t1 += 3600.*24.
                continue

            # Frequency axis
            f = daynoise.f
            nwins.append(np.sum(daynoise.goodwins))

            # Power spectra
            c11_all.append(daynoise.power.c11)
            c22_all.append(daynoise.power.c22)
            cZZ_all.append(daynoise.power.cZZ)
            cPP_all.append(daynoise.power.cPP)

            # Cross spectra
            c12_all.append(daynoise.cross.c12)
            c1Z_all.append(daynoise.cross.c1Z)
            c1P_all.append(daynoise.cross.c1P)
            c2Z_all.append(daynoise.cross.c2Z)
            c2P_all.append(daynoise.cross.c2P)
            cZP_all.append(daynoise.cross.cZP)

            # Rotated spectra
            cHH_all.append(daynoise.rotation.cHH)
            cHZ_all.append(daynoise.rotation.cHZ)
            cHP_all.append(daynoise.rotation.cHP)
            coh_all.append(daynoise.rotation.coh)
            ph_all.append(daynoise.rotation.ph)

            # Coherence
            coh_12_all.append(utils.smooth(utils.coherence(daynoise.cross.c12, daynoise.power.c11, daynoise.power.c22), 50))
            coh_1Z_all.append(utils.smooth(utils.coherence(daynoise.cross.c1Z, daynoise.power.c11, daynoise.power.cZZ), 50))
            coh_1P_all.append(utils.smooth(utils.coherence(daynoise.cross.c1P, daynoise.power.c11, daynoise.power.cPP), 50))
            coh_2Z_all.append(utils.smooth(utils.coherence(daynoise.cross.c2Z, daynoise.power.c22, daynoise.power.cZZ), 50))
            coh_2P_all.append(utils.smooth(utils.coherence(daynoise.cross.c2P, daynoise.power.c22, daynoise.power.cPP), 50))
            coh_ZP_all.append(utils.smooth(utils.coherence(daynoise.cross.cZP, daynoise.power.cZZ, daynoise.power.cPP), 50))

            # Phase
            ph_12_all.append(180./np.pi*utils.phase(daynoise.cross.c12))
            ph_1Z_all.append(180./np.pi*utils.phase(daynoise.cross.c1Z))
            ph_1P_all.append(180./np.pi*utils.phase(daynoise.cross.c1P))
            ph_2Z_all.append(180./np.pi*utils.phase(daynoise.cross.c2Z))
            ph_2P_all.append(180./np.pi*utils.phase(daynoise.cross.c2P))
            ph_ZP_all.append(180./np.pi*utils.phase(daynoise.cross.cZP))

            # Admittance
            ad_12_all.append(utils.smooth(utils.admittance(daynoise.cross.c12, daynoise.power.c11), 50))
            ad_1Z_all.append(utils.smooth(utils.admittance(daynoise.cross.c1Z, daynoise.power.c11), 50))
            ad_1P_all.append(utils.smooth(utils.admittance(daynoise.cross.c1P, daynoise.power.c11), 50))
            ad_2Z_all.append(utils.smooth(utils.admittance(daynoise.cross.c2Z, daynoise.power.c22), 50))
            ad_2P_all.append(utils.smooth(utils.admittance(daynoise.cross.c2P, daynoise.power.c22), 50))
            ad_ZP_all.append(utils.smooth(utils.admittance(daynoise.cross.cZP, daynoise.power.cZZ), 50))

            t1 += 3600.*24.

        # Convert to numpy arrays
        c11_all = np.array(c11_all); c22_all = np.array(c22_all); cZZ_all = np.array(cZZ_all); cPP_all = np.array(cPP_all)
        c12_all = np.array(c12_all); c1Z_all = np.array(c1Z_all); c1P_all = np.array(c1P_all); c2Z_all = np.array(c2Z_all)
        c2P_all = np.array(c2P_all); cZP_all = np.array(cZP_all)
        cHH_all = np.array(cHH_all); cHZ_all = np.array(cHZ_all); cHP_all = np.array(cHP_all)
        coh_all = np.array(coh_all); ph_all = np.array(ph_all)
        coh_12_all = np.array(coh_12_all); coh_1Z_all = np.array(coh_1Z_all); coh_1P_all = np.array(coh_1P_all); coh_2Z_all = np.array(coh_2Z_all)
        coh_2P_all = np.array(coh_2P_all); coh_ZP_all = np.array(coh_ZP_all)
        ph_12_all = np.array(ph_12_all); ph_1Z_all = np.array(ph_1Z_all); ph_1P_all = np.array(ph_1P_all); ph_2Z_all = np.array(ph_2Z_all)
        ph_2P_all = np.array(ph_2P_all); ph_ZP_all = np.array(ph_ZP_all)
        ad_12_all = np.array(ad_12_all); ad_1Z_all = np.array(ad_1Z_all); ad_1P_all = np.array(ad_1P_all); ad_2Z_all = np.array(ad_2Z_all)
        ad_2P_all = np.array(ad_2P_all); ad_ZP_all = np.array(ad_ZP_all)

        nwins = np.array(nwins)

        # Store spectra as objects
        power = Power(c11_all, c22_all, cZZ_all, cPP_all)
        cross = Cross(c12_all, c1Z_all, c1P_all, c2Z_all, c2P_all, cZP_all)
        rotation = Rotation(cHH_all, cHZ_all, cHP_all)

        # Initialize StaNoise object
        stanoise = StaNoise(power, cross, rotation, f, nwins, key=stkey)

        # Store transfer functions as objects for plotting
        coh = Cross(coh_12_all, coh_1Z_all, coh_1P_all, coh_2Z_all, coh_2P_all, coh_ZP_all)
        ph = Cross(ph_12_all, ph_1Z_all, ph_1P_all, ph_2Z_all, ph_2P_all, ph_ZP_all)
        ad = Cross(ad_12_all, ad_1Z_all, ad_1P_all, ad_2Z_all, ad_2P_all, ad_ZP_all)

        # Quality control to identify outliers
        stanoise.QC_sta_spectra(pd=opts.pd, tol=opts.tol, alpha=opts.alpha, 
            fig_QC=opts.fig_QC, debug=opts.debug)

        # Average spectra for good days
        stanoise.average_sta_spectra(fig_average=opte.fig_average, debug=opts.debug)

        if opts.fig_av_cross:
            plot.fig_av_cross(stanoise.f, coh, stanoise.gooddays, 'Coherence', key=stkey, lw=0.5)
            plot.fig_av_cross(stanoise.f, ad, stanoise.gooddays, 'Admittance', key=stkey, lw=0.5)
            plot.fig_av_cross(stanoise.f, ph, stanoise.gooddays, 'Phase', key=stkey, marker=',', lw=0)

        # Save to file
        stanoise.save(fileavst)


if __name__ == "__main__":

    # Run main program
    main()

