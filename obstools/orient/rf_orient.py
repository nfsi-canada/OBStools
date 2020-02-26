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

import sys
import numpy as np
from obspy.core import Stream
import matplotlib.pyplot as plt


def decompose(RF_r, RF_t, t1=0., t2=1., plot_f=False, plot_comps=False):
    """
    Function to decompose radial and transverse receiver function 
    streams into back-azimuth harmonics and determine the main 
    orientation ``azim``, obtained by minimizing the H1' component
    between ``t1`` and ``t2``.

    Parameters
    ----------
    RF_r : :class:`~obspy.core.Stream`
        Stream containing the radial component receiver functions
    RF_t : :class:`~obspy.core.Stream`
        Stream containing the transverse component receiver functions
    t1 : float
        Minimum time over which to calculate ``azcorr`` (sec)
    t2 : float
        Maximum time over which to calculate ``azcorr`` (sec)

    Returns
    -------
    azcorr : float
        Direction (azimuth) along which the Ht_0 harmonic component
        is minimized (between ``t1`` and ``t2``)
    RMS : :class:`~numpy.ndarray`
        Root-mean-square misfit used to determine azcorr
    hr0_rot : :class:`~numpy.ndarray`
        Rotated Hr_0 component
    ht0_rot : :class:`~numpy.ndarray`
        Rotated Ht_0 component

    """

    if not isinstance(RF_r, Stream):
        raise(Exception("Input radial component is not a Stream object"))
    if not isinstance(RF_t, Stream):
        raise(Exception("Input transverse component is not a Stream object"))

    # Some integers
    nbin = len(RF_r)
    nn = len(RF_r[0].data)
    dt = RF_r[0].stats.delta
    daz = 0.1
    naz = int(180./daz)
    deg2rad = np.pi/180.

    # Initialize work arrays
    taxis = np.arange(-nn/2, nn/2)*dt
    trange = np.where((taxis>t1) & (taxis<t2))[0]
    print(trange)
    print(taxis[trange])
    nt = len(trange)
    hr0_rot = np.zeros((nt, naz))
    ht0_rot = np.zeros((nt, naz))
    hr0 = np.zeros(nt); hr1 = np.zeros(nt); hr2 = np.zeros(nt) 
    hr3 = np.zeros(nt); hr4 = np.zeros(nt); meanr = np.zeros(nt)
    ht0 = np.zeros(nt); ht1 = np.zeros(nt); ht2 = np.zeros(nt) 
    ht3 = np.zeros(nt); ht4 = np.zeros(nt); meant = np.zeros(nt)

    # Loop over each depth step
    for ii, it in enumerate(trange):

        # Initialize work arrays
        d_r = np.zeros(nbin)
        d_t = np.zeros(nbin)
        G = np.zeros((nbin, 5))

        # Build arrays and matrices
        for itrace in range(nbin):

            baz = RF_r[itrace].stats.baz
            d_r[itrace] = RF_r[itrace].data[it]
            d_t[itrace] = RF_t[itrace].data[it]
            G[itrace, 0] = 1.0
            G[itrace, 1] = np.cos(deg2rad*baz)
            G[itrace, 2] = np.sin(deg2rad*baz)
            G[itrace, 3] = np.cos(2.*deg2rad*baz)
            G[itrace, 4] = np.sin(2.*deg2rad*baz)

        # Solve using damped least squares
        lam=1.e-25
        m_r = np.linalg.solve(np.dot(G.T, G)+lam*np.identity(G.shape[1]),
            np.dot(G.T, d_r))
        m_t = np.linalg.solve(np.dot(G.T, G)+lam*np.identity(G.shape[1]),
            np.dot(G.T, d_t))

        meanr[ii] = np.mean(d_r)
        hr0[ii] = m_r[0]
        hr1[ii] = m_r[1]
        hr2[ii] = m_r[2]
        hr3[ii] = m_r[3]
        hr4[ii] = m_r[4]

        meant[ii] = np.mean(d_t)
        ht0[ii] = m_t[0]
        ht1[ii] = m_t[1]
        ht2[ii] = m_t[2]
        ht3[ii] = m_t[3]
        ht4[ii] = m_t[4]

        for iaz in range(naz):
            phi = iaz*daz*deg2rad
            hr0_rot[ii, iaz] = np.cos(phi)*m_r[0] + np.sin(phi)*m_t[0]
            ht0_rot[ii, iaz] = -np.sin(phi)*m_r[0] + np.cos(phi)*m_t[0]

    # Minimize misfit of rotated transverse component over specific
    # time range to find azim
    RMS = np.zeros(naz)
    for iaz in range(naz):
        RMS[iaz] = np.sqrt(np.mean(np.square(ht0_rot[:, iaz])))
        # RMS[iaz] = np.sqrt(np.mean(np.square(ht0_rot[indmin:indmax, iaz])))

    # Azimuth of H1
    indaz = np.argmin(RMS)
    azcorr = indaz*daz

    # Resolve ambiguity based on radial component
    if np.mean(hr0_rot[:, indaz]) < 0.:
        azcorr += 180.

    # Rotated components
    phi = deg2rad*azcorr
    meanr_r = np.cos(phi)*meanr + np.sin(phi)*meant
    hr0_r = np.cos(phi)*hr0 + np.sin(phi)*ht0
    hr1_r = np.cos(phi)*hr1 + np.sin(phi)*ht1
    hr2_r = np.cos(phi)*hr2 + np.sin(phi)*ht2
    hr3_r = np.cos(phi)*hr3 + np.sin(phi)*ht3
    hr4_r = np.cos(phi)*hr4 + np.sin(phi)*ht4

    meant_r = -np.sin(phi)*meanr + np.cos(phi)*meant
    ht0_r = -np.sin(phi)*hr0 + np.cos(phi)*ht0
    ht1_r = -np.sin(phi)*hr1 + np.cos(phi)*ht1
    ht2_r = -np.sin(phi)*hr2 + np.cos(phi)*ht2
    ht3_r = -np.sin(phi)*hr3 + np.cos(phi)*ht3
    ht4_r = -np.sin(phi)*hr4 + np.cos(phi)*ht4


    if plot_f:
        fig, ax1 = plt.subplots()
        ax1.set_xlabel("Azimuth (deg)")
        ax1.set_ylabel(r"$f(\phi)$")
        ax1.plot(daz*np.arange(naz), RMS, c='k')
        ax2 = ax1.twinx()
        ax2.plot(daz*np.arange(naz),np.sign(np.mean(hr0_rot,axis=0)),c='tab:blue')
        ax2.tick_params(axis='y', labelcolor='tab:blue')
        ax2.set_ylabel(r"sign(mean($H_{R1}$))", c='tab:blue')
        plt.tight_layout()
        plt.show()

    if plot_comps:
        time = np.linspace(t1, t2, nt)
        fig, ax = plt.subplots(6,4)
        ax[0][0].fill_between(time, 0., meanr, where=meanr>0., facecolor='k', linewidth=0.2)
        ax[1][0].fill_between(time, 0., hr0, where=hr0>0., facecolor='k', linewidth=0.2)
        ax[2][0].fill_between(time, 0., hr1, where=hr1>0., facecolor='k', linewidth=0.2)
        ax[3][0].fill_between(time, 0., hr2, where=hr2>0., facecolor='k', linewidth=0.2)
        ax[4][0].fill_between(time, 0., hr3, where=hr3>0., facecolor='k', linewidth=0.2)
        ax[5][0].fill_between(time, 0., hr4, where=hr4>0., facecolor='k', linewidth=0.2)
        ax[0][0].fill_between(time, 0., meanr, where=meanr<=0., facecolor='r', linewidth=0.2)
        ax[1][0].fill_between(time, 0., hr0, where=hr0<=0., facecolor='r', linewidth=0.2)
        ax[2][0].fill_between(time, 0., hr1, where=hr1<=0., facecolor='r', linewidth=0.2)
        ax[3][0].fill_between(time, 0., hr2, where=hr2<=0., facecolor='r', linewidth=0.2)
        ax[4][0].fill_between(time, 0., hr3, where=hr3<=0., facecolor='r', linewidth=0.2)
        ax[5][0].fill_between(time, 0., hr4, where=hr4<=0., facecolor='r', linewidth=0.2)

        ax[0][1].fill_between(time, 0., meant, where=meant>0., facecolor='k', linewidth=0.2)
        ax[1][1].fill_between(time, 0., ht0, where=ht0>0., facecolor='k', linewidth=0.2)
        ax[2][1].fill_between(time, 0., ht1, where=ht1>0., facecolor='k', linewidth=0.2)
        ax[3][1].fill_between(time, 0., ht2, where=ht2>0., facecolor='k', linewidth=0.2)
        ax[4][1].fill_between(time, 0., ht3, where=ht3>0., facecolor='k', linewidth=0.2)
        ax[5][1].fill_between(time, 0., ht4, where=ht4>0., facecolor='k', linewidth=0.2)
        ax[0][1].fill_between(time, 0., meant, where=meant<=0., facecolor='r', linewidth=0.2)
        ax[1][1].fill_between(time, 0., ht0, where=ht0<=0., facecolor='r', linewidth=0.2)
        ax[2][1].fill_between(time, 0., ht1, where=ht1<=0., facecolor='r', linewidth=0.2)
        ax[3][1].fill_between(time, 0., ht2, where=ht2<=0., facecolor='r', linewidth=0.2)
        ax[4][1].fill_between(time, 0., ht3, where=ht3<=0., facecolor='r', linewidth=0.2)
        ax[5][1].fill_between(time, 0., ht4, where=ht4<=0., facecolor='r', linewidth=0.2)

        ax[0][2].fill_between(time, 0., meanr_r, where=meanr_r>0., facecolor='k', linewidth=0.2)
        ax[1][2].fill_between(time, 0., hr0_r, where=hr0_r>0., facecolor='k', linewidth=0.2)
        ax[2][2].fill_between(time, 0., hr1_r, where=hr1_r>0., facecolor='k', linewidth=0.2)
        ax[3][2].fill_between(time, 0., hr2_r, where=hr2_r>0., facecolor='k', linewidth=0.2)
        ax[4][2].fill_between(time, 0., hr3_r, where=hr3_r>0., facecolor='k', linewidth=0.2)
        ax[5][2].fill_between(time, 0., hr4_r, where=hr4_r>0., facecolor='k', linewidth=0.2)
        ax[0][2].fill_between(time, 0., meanr_r, where=meanr_r<=0., facecolor='r', linewidth=0.2)
        ax[1][2].fill_between(time, 0., hr0_r, where=hr0_r<=0., facecolor='r', linewidth=0.2)
        ax[2][2].fill_between(time, 0., hr1_r, where=hr1_r<=0., facecolor='r', linewidth=0.2)
        ax[3][2].fill_between(time, 0., hr2_r, where=hr2_r<=0., facecolor='r', linewidth=0.2)
        ax[4][2].fill_between(time, 0., hr3_r, where=hr3_r<=0., facecolor='r', linewidth=0.2)
        ax[5][2].fill_between(time, 0., hr4_r, where=hr4_r<=0., facecolor='r', linewidth=0.2)

        ax[0][3].fill_between(time, 0., meant_r, where=meant_r>0., facecolor='k', linewidth=0.2)
        ax[1][3].fill_between(time, 0., ht0_r, where=ht0_r>0., facecolor='k', linewidth=0.2)
        ax[2][3].fill_between(time, 0., ht1_r, where=ht1_r>0., facecolor='k', linewidth=0.2)
        ax[3][3].fill_between(time, 0., ht2_r, where=ht2_r>0., facecolor='k', linewidth=0.2)
        ax[4][3].fill_between(time, 0., ht3_r, where=ht3_r>0., facecolor='k', linewidth=0.2)
        ax[5][3].fill_between(time, 0., ht4_r, where=ht4_r>0., facecolor='k', linewidth=0.2)
        ax[0][3].fill_between(time, 0., meant_r, where=meant_r<=0., facecolor='r', linewidth=0.2)
        ax[1][3].fill_between(time, 0., ht0_r, where=ht0_r<=0., facecolor='r', linewidth=0.2)
        ax[2][3].fill_between(time, 0., ht1_r, where=ht1_r<=0., facecolor='r', linewidth=0.2)
        ax[3][3].fill_between(time, 0., ht2_r, where=ht2_r<=0., facecolor='r', linewidth=0.2)
        ax[4][3].fill_between(time, 0., ht3_r, where=ht3_r<=0., facecolor='r', linewidth=0.2)
        ax[5][3].fill_between(time, 0., ht4_r, where=ht4_r<=0., facecolor='r', linewidth=0.2)

        # Setting the values for all axes.
        plt.setp(ax[0:2], ylim=(-0.0036, 0.0036))
        plt.setp(ax[2:], ylim=(-0.02, 0.02))
        for axis, mode in zip(ax[:, 0], 
            ["Mean", "$H_{1}$", "$H_{2}$", "$H_{3}$", "$H_{4}$", "$H_{5}$"]):
            axis.set_ylabel(mode, size=10)
        ax[0][0].set_title("Radial")
        ax[0][2].set_title("Radial")
        ax[0][1].set_title("Transverse")
        ax[0][3].set_title("Transverse")
        plt.setp(ax[:,1:4], yticklabels=[], yticks=[])
        plt.setp(ax[0:5], xticklabels=[], xticks=[])
        plt.setp(ax[5], xlabel="Time (s)")
        # plt.suptitle("Station: "+RF_r[0].stats.station)
        ax[0, 1].annotate('Misoriented', (0., 0.), xytext=(-12, 60),
                            textcoords='offset points', xycoords='axes fraction',
                            ha='center', va='bottom', size=12)
        ax[0, 3].annotate('Rotated', (0., 0.), xytext=(-12, 60),
                            textcoords='offset points', xycoords='axes fraction',
                            ha='center', va='bottom', size=12)
        plt.show()


    return azcorr, RMS, hr0_rot, ht0_rot


def get_bootstrap(RF_r, RF_t, t1=0., t2=1., plot_hist=True):

    # Bootstrap
    indices = range(len(RF_r))
    isize = int(len(indices)*0.9)
    nboot = 100
    
    az_boot = []
    for iboot in _progressbar(range(nboot), 'Bootstrap sampling: ', 15):
        cc = np.random.choice(indices, size=isize, replace=False)
        RF_r_tmp = Stream()
        RF_t_tmp = Stream()
        for ic in cc:
            RF_r_tmp.append(RF_r[ic])
            RF_t_tmp.append(RF_t[ic])
        azcorr, *_ = decompose(RF_r_tmp, RF_t_tmp, t1, t2)
        az_boot.append(azcorr)

    az_boot = np.array(az_boot)

    if plot_hist:
        plt.hist(az_boot)
        plt.show()

    maz, daz = az_average(az_boot)
    if maz<0.:
        maz += 360.

    return maz, daz

def az_average(azim):
    x = np.cos(azim*np.pi/180.0)
    y = np.sin(azim*np.pi/180.0)
    angles = np.angle(x+1j*y, deg=True)
    angles[angles<0.] += 360.

    phase = np.mean(angles)
    err_phase = np.std(angles)

    return phase, err_phase

def _progressbar(it, prefix="", size=60, file=sys.stdout):
    """
    Show progress bar while looping in for loop

    """

    count = len(it)

    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" %
                   (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()
