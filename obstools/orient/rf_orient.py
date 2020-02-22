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


def decompose(RF_r, RF_t, t1=0., t2=1., plot=True):
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
    nt = len(RF_r[0].data)
    dt = RF_r[0].stats.delta
    daz = 0.1
    naz = int(180./daz)
    deg2rad = np.pi/180.

    # # Define time range over which to calculate azimuth
    # indmin = int(t1/RF_r[0].stats.delta)
    # indmax = int(t2/RF_r[0].stats.delta)

    # Initialize work arrays
    it1 = int(t1/dt)
    it2 = int(t2/dt)
    nt = it2-it1
    hr0_rot = np.zeros((nt, naz))
    ht0_rot = np.zeros((nt, naz))

    # Loop over each depth step
    for ii, it in enumerate(range(it1, it2)):

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

        # Solve system of equations with truncated SVD
        u, s, v = np.linalg.svd(G)
        s[s < 0.001] = 0.
        m_r = np.linalg.solve(s[:, None] * v, u.T.dot(d_r)[:5])
        m_t = np.linalg.solve(s[:, None] * v, u.T.dot(d_t)[:5])

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

    if plot:
        fig, ax1 = plt.subplots()
        ax1.set_xlabel("Azimuth (deg)")
        ax1.set_ylabel(r"$f(\phi)$")
        ax1.plot(daz*np.arange(naz), RMS, c='k')
        ax2 = ax1.twinx()
        ax2.plot(daz*np.arange(naz),np.sign(np.mean(hr0_rot,axis=0)),c='tab:blue')
        ax2.tick_params(axis='y', labelcolor='tab:blue')
        ax2.set_ylabel(r"sign(mean($H_{R1}$))", c='tab:blue')
        plt.show()


    # Resolve ambiguity based on radial component
    if np.mean(hr0_rot[:, indaz]) < 0.:
        azcorr += 180.

    return azcorr, RMS, hr0_rot, ht0_rot


def get_azcorr(RF_r, RF_t, t1=0., t2=1.):

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
        azcorr, *_ = decompose(RF_r_tmp, RF_t_tmp, t1, t2, plot=False)
        az_boot.append(azcorr)

    az_boot = np.array(az_boot)
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
