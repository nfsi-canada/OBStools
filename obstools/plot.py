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
:mod:`~obstools.plot` contains several functions for plotting the results
of the analysis at various final and intermediate steps.

"""

import numpy as np
from matplotlib import pyplot as plt
from obstools import utils, plot


def fig_QC(f, power, gooddays, key=''):
    """
    Function to plot the Quality-Control step of the analysis. This function is used
    in both the `obs_daily_spectra.py` or `obs_clean_spectra.py` scripts.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    power : :class:`~obstools.classes.Power`
        Container for the Power spectra
    gooddays : List
        List of booleans representing whether a window is good (True) or not (False)
    key : str
        String corresponding to the station key under analysis

    """

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
    """
    Function to plot the averaged spectra (those qualified as 'good' in the 
    QC step). This function is used
    in both the `obs_daily_spectra.py` or `obs_clean_spectra.py` scripts.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    power : :class:`~obstools.classes.Power`
        Container for the Power spectra
    bad : :class:`~obstools.classes.Power`
        Container for the *bad* Power spectra
    gooddays : List
        List of booleans representing whether a window is good (True) or not (False)
    key : str
        String corresponding to the station key under analysis

    """

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
    """
    Function to plot the averaged cross-spectra (those qualified as 'good' in the 
    QC step). This function is used in the `obs_daily_spectra.py` script.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    field : :class:`~obstools.classes.Rotation`
        Container for the Power spectra
    gooddays : List
        List of booleans representing whether a window is good (True) or not (False)
    ftype : str
        Type of plot to be displayed. If ftype is Admittance, plot is loglog. Otherwise semilogx
    key : str
        String corresponding to the station key under analysis
    **kwargs : None
        Keyword arguments to modify plot

    """

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
    """
    Function to plot the coherence and phase between the rotated H and Z components, 
    used to characterize the tilt direction.

    Parameters
    ----------
    coh : :mod:`~numpy.ndarray`
        Coherence between rotated H and Z components
    ph : :mod:`~numpy.ndarray`
        Phase between rotated H and Z components
    direc : :mod:`~numpy.ndarray`
        Directions considered in maximizing coherence between H and Z

    """

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


def fig_TF(f, day_trfs, sta_trfs, key=''):
    """
    Function to plot the transfer functions available.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    day_trfs : Dict
        Dictionary containing the transfer functions for the daily averages
    sta_trfs : Dict
        Dictionary containing the transfer functions for the station averages
    key : str
        String corresponding to the station key under analysis

    """

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


def fig_event_raw(evstream, fmin, fmax):
    """
    Function to plot the raw (although bandpassed) seismograms.

    Parameters
    ----------
    evstream : :class:`~obtsools.classes.EventStream`
        Container for the event stream data
    fmin : float
        Low frequency corner (in Hz)
    fmax : float
        High frequency corner (in Hz)

    """

    import matplotlib as mpl
    evstream.sth.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
    evstream.stp.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
    sr = evstream.sth[0].stats.sampling_rate
    taxis = np.arange(0., 7200., 1./sr)
    
    plt.figure(figsize=(6,6))

    plt.subplot(411)
    plt.plot(taxis, evstream.sth[0].data, 'k', lw=0.5)
    plt.title(evstream.key+' '+evstream.tstamp+': H1', fontdict={'fontsize': 8})
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


def fig_event_corrected(evstream, TF_list):
    """
    Function to plot the corrected vertical component seismograms.

    Parameters
    ----------
    evstream : :class:`~obtsools.classes.EventStream`
        Container for the event stream data
    Tf_list : List
        List of Dictionary elements of transfer functions used 
        for plotting the corrected vertical component.

    """

    import matplotlib as mpl
    evstream.sth.filter('bandpass', freqmin=1./150., freqmax = 1./10., corners=2, zerophase=True)
    evstream.stp.filter('bandpass', freqmin=1./150., freqmax = 1./10., corners=2, zerophase=True)
    sr = evstream.sth[0].stats.sampling_rate
    taxis = np.arange(0., 7200., 1./sr)

    plt.figure(figsize=(8,8))
    
    plt.subplot(611)
    plt.plot(taxis, evstream.sth[2].data, 'lightgray', lw=0.5)
    if TF_list['Z1']:
        plt.plot(taxis, evstream.correct['Z1'], 'k', lw=0.5)
    plt.title(evstream.key+' '+evstream.tstamp+': Z1', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True, scilimits=(-3,3))
    plt.xlim((0., 7200.))

    plt.subplot(612)
    plt.plot(taxis, evstream.sth[2].data, 'lightgray', lw=0.5)
    if TF_list['Z2-1']:
        plt.plot(taxis, evstream.correct['Z2-1'], 'k', lw=0.5)
    plt.title(evstream.tstamp+': Z2-1', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True, scilimits=(-3,3))
    plt.xlim((0., 7200.))

    plt.subplot(613)
    plt.plot(taxis, evstream.sth[2].data, 'lightgray', lw=0.5)
    if TF_list['ZP-21']:
        plt.plot(taxis, evstream.correct['ZP-21'], 'k', lw=0.5)
    plt.title(evstream.tstamp+': ZP-21', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True, scilimits=(-3,3))
    plt.xlim((0., 7200.))

    plt.subplot(614)
    plt.plot(taxis, evstream.sth[2].data, 'lightgray', lw=0.5)
    if TF_list['ZH']:
        plt.plot(taxis, evstream.correct['ZH'], 'k', lw=0.5)
    plt.title(evstream.tstamp+': ZH', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True, scilimits=(-3,3))
    plt.xlim((0., 7200.))

    plt.subplot(615)
    plt.plot(taxis, evstream.sth[2].data, 'lightgray', lw=0.5)
    if TF_list['ZP-H']:
        plt.plot(taxis, evstream.correct['ZP-H'], 'k', lw=0.5)
    plt.title(evstream.tstamp+': ZP-H', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True, scilimits=(-3,3))
    plt.xlim((0., 7200.))

    plt.subplot(616)
    plt.plot(taxis, evstream.sth[2].data, 'lightgray', lw=0.5)
    if TF_list['ZP']:
        plt.plot(taxis, evstream.correct['ZP'], 'k', lw=0.5)
    plt.title(evstream.tstamp+': ZP', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True, scilimits=(-3,3))
    plt.xlim((0., 7200.))

    plt.xlabel('Time since earthquake (sec)')
    plt.tight_layout()
    plt.show()


