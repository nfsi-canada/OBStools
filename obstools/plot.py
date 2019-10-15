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

