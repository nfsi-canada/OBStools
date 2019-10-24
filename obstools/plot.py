'''
MODULE obs_plot.py

Set of plotting functions to be used with the obs package for
tilt and compliance correction of vertical component OBS data.

---
Pascal Audet
pascal.audet@uottawa.ca

Last updated: 17 November 2015

'''

import numpy as np
from matplotlib import pyplot as plt
from obstools import utils, plot


def fig_QC(f, power, gooddays, key=''):

    sl_c11 = power.c11
    sl_c22 = power.c22
    sl_cZZ = power.cZZ
    sl_cPP = power.cPP

    plt.figure(6)
    plt.subplot(4,1,1)
    plt.semilogx(f, sl_c11[:,gooddays], 'k', lw=0.5)
    if np.sum(~gooddays)>0:
        plt.semilogx(f, sl_c11[:,~gooddays], 'r', lw=0.5)
    plt.title('H1 component, Station: '+key, fontdict={'fontsize': 8})
    plt.subplot(4,1,2)
    plt.semilogx(f, sl_c22[:,gooddays], 'k', lw=0.5)
    if np.sum(~gooddays)>0:
        plt.semilogx(f, sl_c22[:,~gooddays], 'r', lw=0.5)
    plt.title('H2 component, Station: '+key, fontdict={'fontsize': 8})
    plt.subplot(4,1,3)
    plt.semilogx(f, sl_cZZ[:,gooddays], 'k', lw=0.5)
    if np.sum(~gooddays)>0:
        plt.semilogx(f, sl_cZZ[:,~gooddays], 'r', lw=0.5)
    plt.title('HZ component, Station: '+key, fontdict={'fontsize': 8})
    plt.subplot(4,1,4)
    plt.semilogx(f, sl_cPP[:,gooddays], 'k', lw=0.5)
    if np.sum(~gooddays)>0:
        plt.semilogx(f, sl_cPP[:,~gooddays], 'r', lw=0.5)
    plt.title('HP component, Station: '+key, fontdict={'fontsize': 8})
    plt.tight_layout()
    plt.show()


def fig_average(f, power, bad, gooddays, key=''):

    c11 = power.c11
    c22 = power.c22
    cZZ = power.cZZ
    cPP = power.cPP
    bc11 = bad.c11
    bc22 = bad.c22
    bcZZ = bad.cZZ
    bcPP = bad.cPP

    plt.figure()
    plt.subplot(411)
    plt.semilogx(f, utils.smooth(np.log(c11), 50), 'k', lw=0.5)
    if np.sum(~gooddays)>0:
        plt.semilogx(f, utils.smooth(np.log(bc11), 50), 'r', lw=0.5)
    plt.title('Average H1, Station: '+key, fontdict={'fontsize': 8})
    plt.subplot(412)
    plt.semilogx(f, utils.smooth(np.log(c22), 50), 'k', lw=0.5)
    if np.sum(~gooddays)>0:
        plt.semilogx(f, utils.smooth(np.log(bc22), 50), 'r', lw=0.5)
    plt.title('Average H2, Station: '+key, fontdict={'fontsize': 8})
    plt.subplot(413)
    plt.semilogx(f, utils.smooth(np.log(cZZ), 50), 'k', lw=0.5)
    if np.sum(~gooddays)>0:
        plt.semilogx(f, utils.smooth(np.log(bcZZ), 50), 'r', lw=0.5)
    plt.title('Average HZ, Station: '+key, fontdict={'fontsize': 8})
    plt.subplot(414)
    plt.semilogx(f, utils.smooth(np.log(cPP), 50), 'k', lw=0.5)
    if np.sum(~gooddays)>0:
        plt.semilogx(f, utils.smooth(np.log(bcPP), 50), 'r', lw=0.5)
    plt.title('Average HP, Station: '+key, fontdict={'fontsize': 8})
    plt.tight_layout()
    plt.show()


def fig_av_cross(f, field, gooddays, ftype, key='', **kwargs):

    # Extact field
    if ftype == 'Admittance':
        pp = plt.loglog
    else:
        pp = plt.semilogx

    field12 = field.c12.T
    field1Z = field.c1Z.T
    field1P = field.c1P.T
    field2Z = field.c2Z.T
    field2P = field.c2P.T
    fieldZP = field.cZP.T

    plt.figure(figsize=(6,8))
    plt.subplot(6,1,1)
    pp(f, field12[:,gooddays], color='gray', **kwargs)
    if np.sum(~gooddays)>0:
        pp(f, field12[:,~gooddays], color='r', **kwargs)
    plt.ylabel(ftype, fontdict={'fontsize': 8})
    plt.title(key+' '+ftype+': 12', fontdict={'fontsize': 8})
    plt.subplot(6,1,2)
    pp(f, field1Z[:,gooddays], color='gray', **kwargs)
    if np.sum(~gooddays)>0:
        pp(f, field1Z[:,~gooddays], color='r', **kwargs)
    plt.ylabel(ftype, fontdict={'fontsize': 8})
    plt.title(key+' '+ftype+': 1Z', fontdict={'fontsize': 8})
    plt.subplot(6,1,3)
    pp(f, field1P[:,gooddays], color='gray', **kwargs)
    if np.sum(~gooddays)>0:
        pp(f, field1P[:,~gooddays], color='r', **kwargs)
    plt.ylabel(ftype, fontdict={'fontsize': 8})
    plt.title(key+' '+ftype+': 1P', fontdict={'fontsize': 8})
    plt.subplot(6,1,4)
    pp(f, field2Z[:,gooddays], color='gray', **kwargs)
    if np.sum(~gooddays)>0:
        pp(f, field2Z[:,~gooddays], color='r', **kwargs)
    plt.ylabel(ftype, fontdict={'fontsize': 8})
    plt.title(key+' '+ftype+': 2Z', fontdict={'fontsize': 8})
    plt.subplot(6,1,5)
    pp(f, field2P[:,gooddays], color='gray', **kwargs)
    if np.sum(~gooddays)>0:
        pp(f, field2P[:,~gooddays], color='r', **kwargs)
    plt.ylabel(ftype, fontdict={'fontsize': 8})
    plt.title(key+' '+ftype+': 2P', fontdict={'fontsize': 8})
    plt.subplot(6,1,6)
    pp(f, fieldZP[:,gooddays], color='gray', **kwargs)
    if np.sum(~gooddays)>0:
        pp(f, fieldZP[:,~gooddays], color='r', **kwargs)
    plt.ylabel(ftype, fontdict={'fontsize': 8})
    plt.title(key+' '+ftype+': ZP', fontdict={'fontsize': 8})
    plt.xlabel('Frequency (Hz)', fontdict={'fontsize': 8})
    plt.tight_layout()
    plt.show()


def fig_coh_ph(coh, ph, direc):

    colors = plt.cm.cividis(np.linspace(0,1,coh.shape[1]))

    print(coh.shape)
    if coh.ndim > 1:
        f, (ax1, ax2) = plt.subplots(1,2)
        for i, (co, p) in enumerate(zip(coh, ph)):
            ax1.plot(direc, co, c=colors[i])
            ax2.plot(direc, p*180./np.pi, c=colors[i])
        ax1.set_ylabel('Coherence')
        ax1.set_ylim((0, 1.))
        ax2.set_ylabel('Phase')
        ax1.set_xlabel('Angle from H1')
        ax2.set_xlabel('Angle from H1')
        plt.tight_layout()
        plt.show()
    else:
        plt.figure()
        plt.subplot(121)
        plt.plot(direc, coh, c=colors[0])
        plt.ylim((0, 1.))
        plt.subplot(122)
        plt.plot(direc, ph*180./np.pi, c=colors[0])
        plt.tight_layout()
        plt.show()


def fig_TF(f, day_trfs, list_day, sta_trfs, list_sta, key=''):

    import matplotlib.ticker as mtick
    plt.figure(figsize=(6,8))
    plt.subplot(611)
    for i in range(len(day_trfs)):
        plt.loglog(f, np.abs(day_trfs[i]['ZP']['TF_ZP']), 'gray', lw=0.5)
    plt.loglog(f, np.abs(sta_trfs['ZP']['TF_ZP']), 'k', lw=0.5)
    plt.ylim(1.e-20, 1.e0)
    plt.xlim(1.e-4, 2.5)
    plt.title(key+' Transfer Function: ZP', fontdict={'fontsize': 8})
    plt.subplot(612)
    for i in range(len(day_trfs)):
        plt.loglog(f, np.abs(day_trfs[i]['Z1']['TF_Z1']), 'gray', lw=0.5)
    plt.loglog(f, np.abs(sta_trfs['Z1']['TF_Z1']), 'k', lw=0.5)
    plt.ylim(1.e-5, 1.e5)
    plt.xlim(1.e-4, 2.5)
    plt.title(key+' Transfer Function: Z1', fontdict={'fontsize': 8})
    plt.subplot(613)
    for i in range(len(day_trfs)):
        plt.loglog(f, np.abs(day_trfs[i]['Z2-1']['TF_Z2-1']), 'gray', lw=0.5)
    plt.loglog(f, np.abs(sta_trfs['Z2-1']['TF_Z2-1']), 'k', lw=0.5)
    plt.ylim(1.e-5, 1.e5)
    plt.xlim(1.e-4, 2.5)
    plt.title(key+' Transfer Function: Z2-1', fontdict={'fontsize': 8})
    plt.subplot(614)
    for i in range(len(day_trfs)):
        plt.loglog(f, np.abs(day_trfs[i]['ZP-21']['TF_ZP-21']), 'gray', lw=0.5)
    plt.loglog(f, np.abs(sta_trfs['ZP-21']['TF_ZP-21']), 'k', lw=0.5)
    plt.ylim(1.e-20, 1.e0)
    plt.xlim(1.e-4, 2.5)
    plt.title(key+' Transfer Function: ZP-21', fontdict={'fontsize': 8})
    plt.subplot(615)
    for i in range(len(day_trfs)):
        plt.loglog(f, np.abs(day_trfs[i]['ZH']['TF_ZH']), 'gray', lw=0.5)
    plt.ylim(1.e-10, 1.e10)
    plt.xlim(1.e-4, 2.5)
    plt.title(key+' Transfer Function: ZH', fontdict={'fontsize': 8})
    plt.subplot(616)
    for i in range(len(day_trfs)):
        plt.loglog(f, np.abs(day_trfs[i]['ZP-H']['TF_ZP-H']), 'gray', lw=0.5)
    plt.ylim(1.e-20, 1.e0)
    plt.xlim(1.e-4, 2.5)
    plt.title(key+' Transfer Function: ZP-H', fontdict={'fontsize': 8})
    plt.xlabel('Frequency (Hz)')
    plt.tight_layout()
    plt.show()



def fig_event_raw(evstream):

    import matplotlib as mpl
    evstream.sth.filter('bandpass', freqmin=1./150., freqmax = 1./10., corners=2, zerophase=True)
    evstream.stp.filter('bandpass', freqmin=1./150., freqmax = 1./10., corners=2, zerophase=True)
    sr = evstream.sth[0].stats.sampling_rate
    taxis = np.arange(0., 7200., 1./sr)
    plt.figure(figsize=(6,6))
    plt.subplot(411)
    plt.plot(taxis, evstream.sth[0].data, 'k', lw=0.5)
    plt.title(evstream.tstamp+': H1', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True, scilimits=(-3,3))
    plt.xlim((0., 7200.))
    plt.subplot(412)
    plt.plot(taxis, evstream.sth[1].data, 'k', lw=0.5)
    plt.xlim((0., 7200.))
    plt.title(evstream.tstamp+': H2', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True, scilimits=(-3,3))
    plt.subplot(413)
    plt.plot(taxis, evstream.sth[2].data, 'k', lw=0.5)
    plt.xlim((0., 7200.))
    plt.title(evstream.tstamp+': Z', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True, scilimits=(-3,3))
    plt.subplot(414)
    plt.plot(taxis, evstream.stp[0].data, 'k', lw=0.5)
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True, scilimits=(-3,3))
    plt.xlim((0., 7200.))
    plt.title(evstream.tstamp+': P', fontdict={'fontsize': 8})
    plt.xlabel('Time since earthquake (sec)')
    plt.tight_layout()
    plt.show()


def obs_plot_2traces(tr1, tr2, f1=0.01, f2=0.05, title=None):

    # Get parameters from traces
    nt = tr1.stats.npts
    dt = tr1.stats.delta

    # Time axis
    time = np.arange(nt)*dt

    tt = tr1.copy()
    tt.filter('bandpass', freqmin=f1, freqmax=f2, zerophase=True)
    plt.subplot(211)
    plt.plot(time, tt.data, c='k')

    tt = tr2.copy()
    tt.filter('bandpass', freqmin=f1, freqmax=f2, zerophase=True)
    plt.subplot(212)
    plt.plot(time, tt.data, c='k')

    #if title:
    plt.suptitle(title+tr1.stats.station)

    plt.tight_layout()
    plt.show()
    

def obs_plot_3traces(tr1, tr2, tr3, f1=0.01, f2=0.05, title=None):

    # Get parameters from traces
    nt = tr1.stats.npts
    dt = tr1.stats.delta
    vmax = max(tr1.data.max(),tr2.data.max(),tr3.data.max())

    # Time axis
    time = np.arange(nt)*dt

    tt = tr1.copy()
    tt.filter('bandpass', freqmin=f1, freqmax=f2, zerophase=True)
    plt.subplot(311)
    plt.plot(time, tt.data, c='k')
    plt.ylim(-vmax,vmax)

    tt = tr2.copy()
    tt.filter('bandpass', freqmin=f1, freqmax=f2, zerophase=True)
    plt.subplot(312)
    plt.plot(time, tt.data, c='k')
    plt.ylim(-vmax,vmax)

    tt = tr3.copy()
    tt.filter('bandpass', freqmin=f1, freqmax=f2, zerophase=True)
    plt.subplot(313)
    plt.plot(time, tt.data, c='k')
    plt.ylim(-vmax,vmax)

    #if title:
    plt.suptitle(title+' '+tr1.stats.station)

    plt.tight_layout()
    plt.show()

def obs_plot_trf_tilt(Ad, Co, Ph, eAd, eCo, ePh, freq):

    plt.subplot(231)
    plt.scatter(freq,Co,s=2)
    plt.xlim((0.005,0.1))
    plt.ylim((0., 1.))
    plt.ylabel('Raw coherence')
    plt.xlabel('Frequency (Hz)')

    plt.subplot(232)
    plt.scatter(freq,np.log10(Ad),s=2)
    #plt.errorbar(freq, Ad, yerr=eAd, fmt='.')
    #plt.yscale('log', nonposx='clip')
    plt.xlim((0.005,0.1))
    plt.ylim((-2.5, -1.4))
    plt.ylabel('log10(Admittance)')
    plt.xlabel('Frequency (Hz)')

    plt.subplot(233)
    #plt.scatter(freq,Ph,s=2)
    plt.errorbar(freq, Ph, yerr=ePh, fmt='.', c='k')
    plt.xlim((0.005,0.1))
    plt.ylim((-np.pi,np.pi))
    plt.ylabel('Phase shift (radians)')
    plt.xlabel('Frequency (Hz)')

    plt.suptitle('Transfer functions for tilt')

    plt.tight_layout()
    plt.show()

def obs_plot_trf_compliance(Ad, Co, Ph, eAd, eCo, ePh, freq):

    plt.subplot(231)
    plt.scatter(freq,Co,s=2)
    plt.xlim((0.005,0.1))
    plt.ylim((0., 1.))
    plt.ylabel('Raw coherence')
    plt.xlabel('Frequency (Hz)')

    plt.subplot(232)
    plt.scatter(freq,np.log10(Ad),s=2)
    #plt.errorbar(freq, Ad, yerr=eAd, fmt='.')
    #plt.yscale('log', nonposx='clip')
    plt.xlim((0.005,0.1))
    plt.ylim((-9., -6.))
    plt.ylabel('log10(Admittance)')
    plt.xlabel('Frequency (Hz)')

    plt.subplot(233)
    #plt.scatter(freq,Ph,s=2)
    plt.errorbar(freq, Ph, yerr=ePh, fmt='.', c='k')
    plt.xlim((0.005,0.1))
    plt.ylim((-np.pi,np.pi))
    plt.ylabel('Phase shift (radians)')
    plt.xlabel('Frequency (Hz)')

    plt.suptitle('Transfer functions for compliance')

    plt.tight_layout()
    plt.show()


def obs_plot_coh_dir(dir, coh, tilt):

    if tilt > 180.:
        dir += 180.

    plt.plot(dir,coh, c='k')
    plt.axvline(x=tilt, ls='--', c='k')
    plt.ylabel('Coherence')
    plt.xlabel('Orientation ($^{\circ}$)')
    plt.show()


def obs_plot_trf_analytic(Ad, Ph, eAd, ePh, aAd, aPh, freq, ind):
    
    plt.subplot(221)
    plt.errorbar(freq[ind], Ad[ind], yerr=eAd[ind], fmt='k.')
    plt.plot(freq[ind],aAd[ind],'r')
    #plt.ylim((0.0044,0.0086))
    plt.ylabel('Admittance')
    plt.xlabel('Frequency (Hz)')

    plt.subplot(222)
    plt.errorbar(freq[ind], Ph[ind], yerr=ePh[ind], fmt='k.')
    plt.plot(freq[ind],aPh[ind],'r')
    plt.ylim((-np.pi/10.,np.pi/10.))
    plt.ylabel('Phase shift (radians)')
    plt.xlabel('Frequency (Hz)')

    plt.tight_layout()
    plt.show()

