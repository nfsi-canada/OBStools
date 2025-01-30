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
"""
:mod:`~obstools.plot` contains several functions for plotting the results
of the analysis at various final and intermediate steps.

"""

import numpy as np
from matplotlib import pyplot as plt
from obstools.atacr import utils
from obspy import Trace


def fig_QC(f, power, gooddays, ncomp, key=''):
    """
    Function to plot the Quality-Control step of the analysis. This function
    is used in both the `atacr_daily_spectra` or `atacr_clean_spectra` scripts.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    power : :class:`~obstools.classes.Power`
        Container for the Power spectra
    gooddays : List
        List of booleans representing whether a window is good (True) or not
        (False)
    ncomp : int
        Number of components used in analysis (can be 2, 3 or 4)
    key : str
        String corresponding to the station key under analysis

    """

    sl_c11 = power.c11
    sl_c22 = power.c22
    sl_cZZ = power.cZZ
    sl_cPP = power.cPP

    if ncomp == 2:
        sls = [sl_cZZ, sl_cPP]
        title = ['HZ component, Station: '+key,
                 'HP component, Station: '+key]
    elif ncomp == 3:
        sls = [sl_c11, sl_c22, sl_cZZ]
        title = ['H1 component, Station: '+key,
                 'H2 component, Station: '+key,
                 'HZ component, Station: '+key]
    else:
        sls = [sl_c11, sl_c22, sl_cZZ, sl_cPP]
        title = ['H1 component, Station: '+key,
                 'H2 component, Station: '+key,
                 'HZ component, Station: '+key,
                 'HP component, Station: '+key]

    # Extract only positive frequencies
    faxis = f > 0

    fig = plt.figure(6)
    for i, sl in enumerate(sls):
        ax = fig.add_subplot(ncomp, 1, i+1)
        ax.semilogx(f[faxis], sl[:, gooddays][faxis], 'k', lw=0.5)
        if np.sum(~gooddays) > 0:
            plt.semilogx(f[faxis], sl[:, ~gooddays][faxis], 'r', lw=0.5)
        ax.set_title(title[i], fontdict={'fontsize': 8})
        if i == len(sls)-1:
            plt.xlabel('Frequency (Hz)', fontdict={'fontsize': 8})
    plt.tight_layout()

    return plt


def fig_average(f, power, bad, gooddays, ncomp, key=''):
    """
    Function to plot the averaged spectra (those qualified as 'good' in the
    QC step). This function is used
    in both the `atacr_daily_spectra` or `atacr_clean_spectra` scripts.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    power : :class:`~obstools.classes.Power`
        Container for the Power spectra
    bad : :class:`~obstools.classes.Power`
        Container for the *bad* Power spectra
    gooddays : List
        List of booleans representing whether a window is good (True) or not
        (False)
    ncomp : int
        Number of components used in analysis (can be 2, 3 or 4)
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

    if ncomp == 2:
        ccs = [cZZ, cPP]
        bcs = [bcZZ, bcPP]
        title = ['Average HZ, Station: '+key,
                 'Average HP, Station: '+key]
    elif ncomp == 3:
        ccs = [c11, c22, cZZ]
        bcs = [bc11, bc22, bcZZ]
        title = ['Average H1, Station: '+key,
                 'Average H2, Station: '+key,
                 'Average HZ, Station: '+key]
    else:
        ccs = [c11, c22, cZZ, cPP]
        bcs = [bc11, bc22, bcZZ, bcPP]
        title = ['Average H1, Station: '+key,
                 'Average H2, Station: '+key,
                 'Average HZ, Station: '+key,
                 'Average HP, Station: '+key]

    # Extract only positive frequencies
    faxis = f > 0

    plt.figure()
    for i, (cc, bc) in enumerate(zip(ccs, bcs)):
        ax = plt.subplot(ncomp, 1, i+1)
        ax.semilogx(
            f[faxis], utils.smooth(np.log(cc)[faxis], 50), 'k', lw=0.5)
        if np.sum(~gooddays) > 0:
            ax.semilogx(
                f[faxis], utils.smooth(np.log(bc)[faxis], 50), 'r', lw=0.5)
        ax.set_title(title[i], fontdict={'fontsize': 8})
        if i == len(ccs)-1:
            plt.xlabel('Frequency (Hz)', fontdict={'fontsize': 8})
    plt.tight_layout()

    return plt


def fig_av_cross(f, field, gooddays, ftype, ncomp, key='',
                 save=False, fname='', form='png', **kwargs):
    """
    Function to plot the averaged cross-spectra (those qualified as 'good' in
    the QC step). This function is used in the `atacr_daily_spectra` script.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    field : :class:`~obstools.classes.Rotation`
        Container for the Power spectra
    gooddays : List
        List of booleans representing whether a window is good (True) or not
        (False)
    ftype : str
        Type of plot to be displayed. If ftype is Admittance, plot is loglog.
        Otherwise semilogx
    key : str
        String corresponding to the station key under analysis
    **kwargs : None
        Keyword arguments to modify plot

    """

    # Extract only positive frequencies
    faxis = f > 0

    if ncomp == 2:
        fieldZP = field.cZP.T
        fields = [fieldZP]
        title = [': ZP']
        fig = plt.figure(figsize=(6, 2.667))
    elif ncomp == 3:
        field12 = field.c12.T
        field1Z = field.c1Z.T
        field2Z = field.c2Z.T
        fields = [field12, field1Z, field2Z]
        title = [': 12', ': 1Z', ': 2Z']
        fig = plt.figure(figsize=(6, 4))
    else:
        fieldZP = field.cZP.T
        field12 = field.c12.T
        field1Z = field.c1Z.T
        field2Z = field.c2Z.T
        field1P = field.c1P.T
        field2P = field.c2P.T
        fields = [field12, field1Z, field1P, field2Z, field2P, fieldZP]
        title = [': 12', ': 1Z', ': 1P', ': 2Z', ': 2P', ': ZP']
        fig = plt.figure(figsize=(6, 8))

    for i, field in enumerate(fields):
        ax = fig.add_subplot(len(fields), 1, i+1)
        # Extact field
        if ftype == 'Admittance':
            ax.loglog(
                f[faxis], field[:, gooddays][faxis], color='gray', **kwargs)
            if np.sum(~gooddays) > 0:
                ax.loglog(
                    f[faxis], field[:, ~gooddays][faxis], color='r', **kwargs)
        else:
            ax.semilogx(
                f[faxis], field[:, gooddays][faxis], color='gray', **kwargs)
            if np.sum(~gooddays) > 0:
                ax.semilogx(
                    f[faxis], field[:, ~gooddays][faxis], color='r', **kwargs)
        plt.ylabel(ftype, fontdict={'fontsize': 8})
        plt.title(key+' '+ftype+title[i], fontdict={'fontsize': 8})
        if i == len(fields)-1:
            plt.xlabel('Frequency (Hz)', fontdict={'fontsize': 8})

    plt.tight_layout()

    return plt


def fig_coh_ph(coh, ph, direc):
    """
    Function to plot the coherence and phase between the rotated H and Z
    components, used to characterize the tilt direction.

    Parameters
    ----------
    coh : :mod:`~numpy.ndarray`
        Coherence between rotated H and Z components
    ph : :mod:`~numpy.ndarray`
        Phase between rotated H and Z components
    direc : :mod:`~numpy.ndarray`
        Directions considered in maximizing coherence between H and Z

    """

    colors = plt.cm.cividis(np.linspace(0, 1, coh.shape[0]))

    if coh.ndim > 1:
        f, (ax1, ax2) = plt.subplots(1, 2)
        for i, (co, p) in enumerate(zip(coh, ph)):
            ax1.plot(direc, co, c=colors[i])
            ax2.plot(direc, p*180./np.pi, c=colors[i])
        ax1.set_ylabel('Coherence')
        ax1.set_ylim((0, 1.))
        ax2.set_ylabel('Phase')
        ax1.set_xlabel('Angle from H1')
        ax2.set_xlabel('Angle from H1')
        plt.tight_layout()

    else:
        plt.figure()
        plt.subplot(121)
        plt.plot(direc, coh, c=colors[0])
        plt.ylim((0, 1.))
        plt.subplot(122)
        plt.plot(direc, ph*180./np.pi, c=colors[0])
        plt.tight_layout()

    return plt


def fig_TF(f, day_trfs, day_list, sta_trfs, sta_list, skey=''):
    """
    Function to plot the transfer functions available.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    day_trfs : Dict
        Dictionary containing the transfer functions for the daily averages
    day_list : Dict
        Dictionary containing the list of daily transfer functions
    sta_trfs : Dict
        Dictionary containing the transfer functions for the station averages
    sta_list : Dict
        Dictionary containing the list of average transfer functions
    skey : str
        String corresponding to the station key under analysis

    """

    import matplotlib.ticker as mtick

    # Extract only positive frequencies
    faxis = f > 0

    # Get max number of TFs to plot
    ntf = max(sum(day_list.values()), sum(sta_list.values()))

    # Define all possible compbinations
    tf_list = {'ZP': True, 'Z1': True, 'Z2-1': True,
               'ZP-21': True, 'ZH': True, 'ZP-H': True}

    if ntf == 1:
        fig = plt.figure(figsize=(6, 1.75))
    else:
        fig = plt.figure(figsize=(6, 1.33333333*ntf))

    j = 1
    for key in tf_list:

        if not day_list[key] and not sta_list[key]:
            continue

        ax = fig.add_subplot(ntf, 1, j)

        if day_list[key]:
            for i in range(len(day_trfs)):
                ax.loglog(
                    f[faxis],
                    np.abs(day_trfs[i][key]['TF_'+key][faxis]),
                    'gray', lw=0.5)
        if sta_list[key]:
            ax.loglog(
                f[faxis],
                np.abs(sta_trfs[key]['TF_'+key][faxis]),
                'k', lw=0.5)
        if key == 'ZP':
            ax.set_ylim(1.e-20, 1.e0)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(skey+' Transfer Function: ZP',
                         fontdict={'fontsize': 8})
        elif key == 'Z1':
            ax.set_ylim(1.e-5, 1.e5)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(skey+' Transfer Function: Z1',
                         fontdict={'fontsize': 8})
        elif key == 'Z2-1':
            ax.set_ylim(1.e-5, 1.e5)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(skey+' Transfer Function: Z2-1',
                         fontdict={'fontsize': 8})
        elif key == 'ZP-21':
            ax.set_ylim(1.e-20, 1.e0)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(skey+' Transfer Function: ZP-21',
                         fontdict={'fontsize': 8})
        elif key == 'ZH':
            ax.set_ylim(1.e-10, 1.e10)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(skey+' Transfer Function: ZH',
                         fontdict={'fontsize': 8})
        elif key == 'ZP-H':
            ax.set_ylim(1.e-20, 1.e0)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(skey+' Transfer Function: ZP-H',
                         fontdict={'fontsize': 8})

        j += 1

    ax.set_xlabel('Frequency (Hz)')
    plt.tight_layout()

    return plt


def fig_comply(f, day_comps, day_list, sta_comps, sta_list, skey=None,
    elev=-1000., f_0=None):
    """
    Function to plot the transfer functions available.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    day_comps : Dict
        Dictionary containing the compliance functions for the daily averages
    day_list : Dict
        Dictionary containing the list of daily transfer functions
    sta_comps : Dict
        Dictionary containing the compliance functions for the station averages
    sta_list : Dict
        Dictionary containing the list of average transfer functions
    skey : str
        String corresponding to the station key under analysis
    elev : float
        Station elevation in meters (OBS stations have negative elevations)
    f_0 : float
        Lowest frequency to consider in plot (Hz)

    """

    import matplotlib.ticker as mtick
    import matplotlib.pyplot as plt

    # Extract only positive frequencies
    faxis = f > 0

    # Positive station elevation for frequency limit calc
    elev = -1.*elev

    # Calculate theoretical frequency limit for infra-gravity waves
    f_c = np.sqrt(9.81/np.pi/elev)/2.

    # Define all possible combinations
    comp_list = {'ZP': True, 'ZP-21': True, 'ZP-H': True}

    # Get max number of subplot
    nkeys_day = sum(day_list[key] for key in comp_list)
    nkeys_sta = sum(sta_list[key] for key in comp_list)
    ncomps = max(nkeys_day, nkeys_sta)

    if ncomps == 1:
        fig = plt.figure(figsize=(6, 1.75))
    else:
        fig = plt.figure(figsize=(6, 1.33333333*ncomps))

    for j, key in enumerate(comp_list):

        if not day_list[key] and not sta_list[key]:
            continue

        ax = fig.add_subplot(ncomps, 2, j*2+1)
        ax.tick_params(labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)

        if day_list[key]:
            compliance_list = []
            coherence_list = []
            for i in range(len(day_comps)):
                compliance = np.abs(day_comps[i][key][0])
                coherence = np.abs(day_comps[i][key][1])
                if not np.isnan(compliance).any():
                    compliance_list.append(compliance)
                    coherence_list.append(coherence)
            compliance_mean = np.mean(np.array(compliance_list), axis=0)
            compliance_std = np.std(np.array(compliance_list), axis=0)
            coherence_mean = np.mean(np.array(coherence_list), axis=0)
            coherence_std = np.std(np.array(coherence_list), axis=0)

            ax.fill_between(
                f[faxis], 
                compliance_mean[faxis]-compliance_std[faxis], 
                compliance_mean[faxis]+compliance_std[faxis], 
                fc='royalblue', alpha=0.3, label=r'$\pm$ Std daily'
                )
            ax.plot(
                f[faxis], compliance_mean[faxis], c='royalblue', 
                lw=0.5, label='Mean daily')
            ax.set_xlim(f_0, f_c)
            ytop = 1.2*np.max(compliance_mean[(f > f_0) & (f < f_c)])
            ybot = 0/8*np.min(compliance_mean[(f > f_0) & (f < f_c)])
            ax.set_ylim(ybot, ytop)

        if sta_list[key]:
            for i in range(len(sta_comps)):
                compliance = np.abs(sta_comps[i][key][0])
                ax.plot(
                    f[faxis],
                    compliance[faxis],
                    'red', lw=0.5, alpha=0.5,
                    label='Sta average')

        if key == 'ZP':
            ax.set_title(skey+' Compliance: ZP',
                         fontdict={'fontsize': 8})
        elif key == 'ZP-21':
            ax.set_title(skey+' Compliance: ZP-21',
                         fontdict={'fontsize': 8})
        elif key == 'ZP-H':
            ax.set_title(skey+' Compliance: ZP-H',
                         fontdict={'fontsize': 8})

        if f_0:
            ax.axvline(f_0, ls='--', c='k', lw=0.75)
        ax.axvline(f_c, ls='--', c='k', lw=0.75)

        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), fontsize=6)

        ax = fig.add_subplot(ncomps, 2, j*2+2)
        ax.tick_params(labelsize=8)

        if day_list[key]:
            # for i in range(len(day_comps)):
            ax.fill_between(
                f[faxis],
                coherence_mean[faxis]-coherence_std[faxis],
                coherence_mean[faxis]+coherence_std[faxis],
                fc='royalblue', alpha=0.3
                )
            ax.plot(
                f[faxis],
                coherence_mean[faxis],
                c='royalblue', lw=0.75)
        if sta_list[key]:
            for i in range(len(sta_comps)):
                ax.plot(
                    f[faxis], 
                    np.abs(sta_comps[i][key][1][faxis]),
                    'red', lw=0.5, alpha=0.5)
        ax.set_xscale('log')

        if key == 'ZP':
            ax.set_title(skey+' Coherence: ZP',
                         fontdict={'fontsize': 8})
        elif key == 'ZP-21':
            ax.set_title(skey+' Coherence: ZP-21',
                         fontdict={'fontsize': 8})
        elif key == 'ZP-H':
            ax.set_title(skey+' Coherence: ZP-H',
                         fontdict={'fontsize': 8})

        if f_0:
            ax.axvline(f_0, ls='--', c='k', lw=0.75)
        ax.axvline(f_c, ls='--', c='k', lw=0.75)

    axes = plt.gcf().get_axes()
    axes[-2].set_xlabel('Frequency (Hz)', fontsize=8)
    axes[-1].set_xlabel('Frequency (Hz)', fontsize=8)

    plt.tight_layout()

    return plt


def fig_event_raw(evstream, fmin=1./150., fmax=2.):
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

    from obspy import Stream

    # Unpack traces
    tr1 = evstream.tr1.copy()
    tr2 = evstream.tr2.copy()
    trZ = evstream.trZ.copy()
    trP = evstream.trP.copy()
    st = Stream(traces=[tr for tr in [tr1, tr2, trZ, trP] if np.any(tr.data)])
    st.filter(
        'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
    sr = trZ.stats.sampling_rate

    taxis = np.arange(0., trZ.stats.npts/sr, 1./sr)

    fig = plt.figure(figsize=(6, 6))

    ax = fig.add_subplot(4, 1, 1)
    ax.plot(taxis, trZ.data, 'k', lw=0.5)
    ax.set_title(evstream.key + ' ' + evstream.tstamp +
                 ': Z', fontdict={'fontsize': 8})
    ax.ticklabel_format(axis='y', style='sci', useOffset=True,
                        scilimits=(-3, 3))
    ax.set_xlim((0., trZ.stats.npts/sr))

    if len(tr1.data > 0):
        ax = fig.add_subplot(4, 1, 2)
        ax.plot(taxis, tr1.data, 'k', lw=0.5)
        ax.set_xlim((0., 7200.))
        ax.set_title(evstream.tstamp + ': 1', fontdict={'fontsize': 8})
        ax.ticklabel_format(axis='y', style='sci', useOffset=True,
                            scilimits=(-3, 3))

        ax = fig.add_subplot(4, 1, 3)
        ax.plot(taxis, tr2.data, 'k', lw=0.5)
        ax.set_xlim((0., trZ.stats.npts/sr))
        ax.set_title(evstream.tstamp + ': 2', fontdict={'fontsize': 8})
        ax.ticklabel_format(axis='y', style='sci', useOffset=True,
                            scilimits=(-3, 3))

    if len(trP.data > 0):
        if len(tr1.data > 0):
            ax = fig.add_subplot(4, 1, 4)
        else:
            ax = fig.add_subplot(4, 1, 2)
        ax.plot(taxis, trP.data, 'k', lw=0.5)
        ax.ticklabel_format(axis='y', style='sci', useOffset=True,
                            scilimits=(-3, 3))
        ax.set_xlim((0., trZ.stats.npts/sr))
        ax.set_title(evstream.tstamp + ': P', fontdict={'fontsize': 8})

    plt.xlabel('Time since earthquake (sec)')
    plt.tight_layout()

    return plt


def fig_event_corrected(evstream, TF_list, fmin=1./150., fmax=2.):
    """
    Function to plot the corrected vertical component seismograms.

    Parameters
    ----------
    evstream : :class:`~obtsools.classes.EventStream`
        Container for the event stream data
    Tf_list : list
        List of Dictionary elements of transfer functions used
        for plotting the corrected vertical component.

    """

    # Unpack vertical trace and filter
    trZ = evstream.trZ.copy()
    trZ.filter(
        'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
    sr = trZ.stats.sampling_rate
    taxis = np.arange(0., trZ.stats.npts/sr, 1./sr)

    plt.figure(figsize=(8, 8))

    plt.subplot(611)
    plt.plot(
        taxis, trZ.data, 'lightgray', lw=0.5)
    if TF_list['Z1']:
        tr = Trace(
            data=evstream.correct['Z1'],
            header=trZ.stats).filter(
            'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
        plt.plot(taxis, tr.data, 'k', lw=0.5)
    plt.title(evstream.key + ' ' + evstream.tstamp +
              ': Z1', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                               scilimits=(-3, 3))
    plt.xlim((0., trZ.stats.npts/sr))

    plt.subplot(612)
    plt.plot(
        taxis, trZ.data, 'lightgray', lw=0.5)
    if TF_list['Z2-1']:
        tr = Trace(
            data=evstream.correct['Z2-1'],
            header=trZ.stats).filter(
            'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
        plt.plot(taxis, tr.data, 'k', lw=0.5)
    plt.title(evstream.tstamp + ': Z2-1', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                               scilimits=(-3, 3))
    plt.xlim((0., trZ.stats.npts/sr))

    plt.subplot(613)
    plt.plot(
        taxis, trZ.data, 'lightgray', lw=0.5)
    if TF_list['ZP-21']:
        tr = Trace(
            data=evstream.correct['ZP-21'],
            header=trZ.stats).filter(
            'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
        plt.plot(taxis, tr.data, 'k', lw=0.5)
    plt.title(evstream.tstamp + ': ZP-21', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                               scilimits=(-3, 3))
    plt.xlim((0., trZ.stats.npts/sr))

    plt.subplot(614)
    plt.plot(
        taxis, trZ.data, 'lightgray', lw=0.5)
    if TF_list['ZH']:
        tr = Trace(
            data=evstream.correct['ZH'],
            header=trZ.stats).filter(
            'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
        plt.plot(taxis, tr.data, 'k', lw=0.5)
    plt.title(evstream.tstamp + ': ZH', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                               scilimits=(-3, 3))
    plt.xlim((0., trZ.stats.npts/sr))

    plt.subplot(615)
    plt.plot(
        taxis, trZ.data, 'lightgray', lw=0.5)
    if TF_list['ZP-H']:
        tr = Trace(
            data=evstream.correct['ZP-H'],
            header=trZ.stats).filter(
            'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
        plt.plot(taxis, tr.data, 'k', lw=0.5)
    plt.title(evstream.tstamp + ': ZP-H', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                               scilimits=(-3, 3))
    plt.xlim((0., trZ.stats.npts/sr))

    plt.subplot(616)
    plt.plot(
        taxis, trZ.data, 'lightgray', lw=0.5)
    if TF_list['ZP']:
        tr = Trace(
            data=evstream.correct['ZP'],
            header=trZ.stats).filter(
            'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
        plt.plot(taxis, tr.data, 'k', lw=0.5)
    plt.title(evstream.tstamp + ': ZP', fontdict={'fontsize': 8})
    plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                               scilimits=(-3, 3))
    plt.xlim((0., trZ.stats.npts/sr))

    plt.xlabel('Time since earthquake (sec)')
    plt.tight_layout()

    return plt
