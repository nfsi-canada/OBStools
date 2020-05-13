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
import os
import sys
import numpy as np
import pickle
import stdb
from obstools.atacr import StaNoise, Power, Cross, Rotation
from obstools.atacr import utils, arguments, plot


def main():

    # Run Input Parser
    args = arguments.get_cleanspec_arguments()

    # Load Database
    db = stdb.io.load_db(fname=args.indb)

    # Construct station key loop
    allkeys = db.keys()
    sorted(allkeys)

    # Extract key subset
    if len(args.stkeys) > 0:
        stkeys = []
        for skey in args.stkeys:
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
        if not os.path.isdir(specpath):
            print("Path to "+specpath+" doesn`t exist - aborting")
            sys.exit()

        # Path where average spectra will be saved
        avstpath = 'AVG_STA/' + stkey + '/'
        if not os.path.isdir(avstpath):
            print("Path to "+avstpath+" doesn`t exist - creating it")
            os.makedirs(avstpath)

        # Path where plots will be saved
        if args.saveplot:
            plotpath = avstpath + 'PLOTS/'
            if not os.path.isdir(plotpath):
                os.makedirs(plotpath)
        else:
            plotpath = False

        # Get catalogue search start time
        if args.startT is None:
            tstart = sta.startdate
        else:
            tstart = args.startT

        # Get catalogue search end time
        if args.endT is None:
            tend = sta.enddate
        else:
            tend = args.endT

        if tstart > sta.enddate or tend < sta.startdate:
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
        print()
        print("|===============================================|")
        print("|===============================================|")
        print("|                   {0:>8s}                    |".format(
            sta.station))
        print("|===============================================|")
        print("|===============================================|")
        print("|  Station: {0:>2s}.{1:5s}                            |".format(
            sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}  |".format(
            sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(
            sta.longitude, sta.latitude))
        print("|      Start time: {0:19s}          |".format(
            sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|      End time:   {0:19s}          |".format(
            sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|-----------------------------------------------|")

        # Filename for output average spectra
        dstart = str(tstart.year).zfill(4)+'.'+str(tstart.julday).zfill(3)+'-'
        dend = str(tend.year).zfill(4)+'.'+str(tend.julday).zfill(3)+'.'
        fileavst = avstpath + dstart + dend + 'avg_sta.pkl'

        if os.path.exists(fileavst):
            if not args.ovr:
                print("*   -> file "+fileavst+" exists - continuing")
                continue

        # Containers for power and cross spectra
        coh_all = []
        ph_all = []
        coh_12_all = []
        coh_1Z_all = []
        coh_1P_all = []
        coh_2Z_all = []
        coh_2P_all = []
        coh_ZP_all = []
        ph_12_all = []
        ph_1Z_all = []
        ph_1P_all = []
        ph_2Z_all = []
        ph_2P_all = []
        ph_ZP_all = []
        ad_12_all = []
        ad_1Z_all = []
        ad_1P_all = []
        ad_2Z_all = []
        ad_2P_all = []
        ad_ZP_all = []
        nwins = []

        t1 = tstart

        # Initialize StaNoise object
        stanoise = StaNoise()

        # Loop through each day withing time range
        while t1 < tend:

            year = str(t1.year).zfill(4)
            jday = str(t1.julday).zfill(3)

            tstamp = year+'.'+jday+'.'
            filespec = specpath + tstamp + 'spectra.pkl'

            # Load file if it exists
            if os.path.exists(filespec):
                print()
                print(
                    "*******************************************" +
                    "*****************")
                print('* Calculating noise spectra for key ' +
                      stkey+' and day '+year+'.'+jday)
                print("*   -> file "+filespec+" found - loading")
                file = open(filespec, 'rb')
                daynoise = pickle.load(file)
                file.close()
                stanoise += daynoise
            else:
                t1 += 3600.*24.
                continue

            coh_all.append(daynoise.rotation.coh)
            ph_all.append(daynoise.rotation.ph)

            # Coherence
            coh_12_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.c12,
                        daynoise.power.c11,
                        daynoise.power.c22), 50))
            coh_1Z_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.c1Z,
                        daynoise.power.c11,
                        daynoise.power.cZZ), 50))
            coh_1P_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.c1P,
                        daynoise.power.c11,
                        daynoise.power.cPP), 50))
            coh_2Z_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.c2Z,
                        daynoise.power.c22,
                        daynoise.power.cZZ), 50))
            coh_2P_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.c2P,
                        daynoise.power.c22,
                        daynoise.power.cPP), 50))
            coh_ZP_all.append(
                utils.smooth(
                    utils.coherence(
                        daynoise.cross.cZP,
                        daynoise.power.cZZ,
                        daynoise.power.cPP), 50))

            # Phase
            try:
                ph_12_all.append(
                    180./np.pi*utils.phase(daynoise.cross.c12))
            except:
                ph_12_all.append(None)
            try:
                ph_1Z_all.append(
                    180./np.pi*utils.phase(daynoise.cross.c1Z))
            except:
                ph_1Z_all.append(None)
            try:
                ph_1P_all.append(
                    180./np.pi*utils.phase(daynoise.cross.c1P))
            except:
                ph_1P_all.append(None)
            try:
                ph_2Z_all.append(
                    180./np.pi*utils.phase(daynoise.cross.c2Z))
            except:
                ph_2Z_all.append(None)
            try:
                ph_2P_all.append(
                    180./np.pi*utils.phase(daynoise.cross.c2P))
            except:
                ph_2P_all.append(None)
            try:
                ph_ZP_all.append(
                    180./np.pi*utils.phase(daynoise.cross.cZP))
            except:
                ph_ZP_all.append(None)

            # Admittance
            ad_12_all.append(utils.smooth(utils.admittance(
                daynoise.cross.c12, daynoise.power.c11), 50))
            ad_1Z_all.append(utils.smooth(utils.admittance(
                daynoise.cross.c1Z, daynoise.power.c11), 50))
            ad_1P_all.append(utils.smooth(utils.admittance(
                daynoise.cross.c1P, daynoise.power.c11), 50))
            ad_2Z_all.append(utils.smooth(utils.admittance(
                daynoise.cross.c2Z, daynoise.power.c22), 50))
            ad_2P_all.append(utils.smooth(utils.admittance(
                daynoise.cross.c2P, daynoise.power.c22), 50))
            ad_ZP_all.append(utils.smooth(utils.admittance(
                daynoise.cross.cZP, daynoise.power.cZZ), 50))

            t1 += 3600.*24.

        # Convert to numpy arrays
        coh_all = np.array(coh_all)
        ph_all = np.array(ph_all)
        coh_12_all = np.array(coh_12_all)
        coh_1Z_all = np.array(coh_1Z_all)
        coh_1P_all = np.array(coh_1P_all)
        coh_2Z_all = np.array(coh_2Z_all)
        coh_2P_all = np.array(coh_2P_all)
        coh_ZP_all = np.array(coh_ZP_all)
        ph_12_all = np.array(ph_12_all)
        ph_1Z_all = np.array(ph_1Z_all)
        ph_1P_all = np.array(ph_1P_all)
        ph_2Z_all = np.array(ph_2Z_all)
        ph_2P_all = np.array(ph_2P_all)
        ph_ZP_all = np.array(ph_ZP_all)
        ad_12_all = np.array(ad_12_all)
        ad_1Z_all = np.array(ad_1Z_all)
        ad_1P_all = np.array(ad_1P_all)
        ad_2Z_all = np.array(ad_2Z_all)
        ad_2P_all = np.array(ad_2P_all)
        ad_ZP_all = np.array(ad_ZP_all)

        # Store transfer functions as objects for plotting
        coh = Cross(coh_12_all, coh_1Z_all, coh_1P_all,
                    coh_2Z_all, coh_2P_all, coh_ZP_all)
        ph = Cross(ph_12_all, ph_1Z_all, ph_1P_all,
                   ph_2Z_all, ph_2P_all, ph_ZP_all)
        ad = Cross(ad_12_all, ad_1Z_all, ad_1P_all,
                   ad_2Z_all, ad_2P_all, ad_ZP_all)

        # Quality control to identify outliers
        stanoise.QC_sta_spectra(pd=args.pd, tol=args.tol, alpha=args.alpha,
                                fig_QC=args.fig_QC, debug=args.debug,
                                save=plotpath, form=args.form)

        # Average spectra for good days
        stanoise.average_sta_spectra(
            fig_average=args.fig_average,
            save=plotpath, form=args.form)

        if args.fig_av_cross:
            fname = stkey + '.' + 'av_cross'
            plot.fig_av_cross(stanoise.f, coh, stanoise.gooddays,
                              'Coherence', stanoise.ncomp, key=stkey, lw=0.5,
                              save=plotpath, fname=fname, form=args.form)
            plot.fig_av_cross(stanoise.f, ad, stanoise.gooddays,
                              'Admittance', stanoise.ncomp, key=stkey, lw=0.5,
                              save=plotpath, fname=fname, form=args.form)
            plot.fig_av_cross(stanoise.f, ph, stanoise.gooddays,
                              'Phase', stanoise.ncomp, key=stkey, marker=',', lw=0,
                              save=plotpath, fname=fname, form=args.form)

        if args.fig_coh_ph and stanoise.direc.any():
            fname = stkey + '.' + 'coh_ph'
            plot.fig_coh_ph(coh_all, ph_all, stanoise.direc,
                save=plotpath, fname=fname, form=args.form)

        # Save to file
        stanoise.save(fileavst)


if __name__ == "__main__":

    # Run main program
    main()
