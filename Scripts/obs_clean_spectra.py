#!/usr/bin/env python

# Import modules and functions
import os
import numpy as np
from obspy import UTCDateTime
import pickle
from obstools import StaNoise, Power, Cross, Rotation
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

    # Path where average spectra will be saved
    avstpath = 'AVG_STA/' + sta_key[0] + '/'
    if not os.path.isdir(avstpath): 
        print('Path to '+avstpath+' doesn`t exist - creating it')
        os.makedirs(avstpath)

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
    # Lopp through each day withing time range
    while t1 < tend:

        year = str(t1.year).zfill(4)
        jday = str(t1.julday).zfill(3)

        print('Calculating noise spectra for key '+sta_key[0]+' and day '+year+'.'+jday)
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
    stanoise = StaNoise(power, cross, rotation, f, nwins, key=sta_key[0])

    # Store transfer functions as objects for plotting
    coh = Cross(coh_12_all, coh_1Z_all, coh_1P_all, coh_2Z_all, coh_2P_all, coh_ZP_all)
    ph = Cross(ph_12_all, ph_1Z_all, ph_1P_all, ph_2Z_all, ph_2P_all, ph_ZP_all)
    ad = Cross(ad_12_all, ad_1Z_all, ad_1P_all, ad_2Z_all, ad_2P_all, ad_ZP_all)

    # Quality control to identify outliers
    stanoise.QC_sta_spectra(pd=[0.004, 2.0], fig_QC=False, debug=False)

    # Average spectra for good days
    stanoise.average_sta_spectra(fig_average=False, debug=False)

    # plot.fig_coh_ph(coh_all, ph_all, np.arange(0., 360., 10.))

    # if fig_av_cross:
    plot.fig_av_cross(stanoise.f, coh, stanoise.gooddays, 'Coherence', key=sta_key[0], lw=0.5)
    plot.fig_av_cross(stanoise.f, ad, stanoise.gooddays, 'Admittance', key=sta_key[0], lw=0.5)
    plot.fig_av_cross(stanoise.f, ph, stanoise.gooddays, 'Phase', key=sta_key[0], marker=',', lw=0)

    # Save to file
    stanoise.save(fileavst)


if __name__ == "__main__":

    # Run main program
    main()

