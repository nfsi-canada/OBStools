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

import numpy as np
from obspy.core import Stream


def find_azcorr(RF_r, RF_t, t1=0., t2=1.):
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

    print()
    print('Decomposing receiver functions into baz harmonics')

    if not isinstance(RF_r, obspy.core.Stream):
        raise(Exception("Input radial component is not a Stream object"))
    if not isinstance(RF_t, obspy.core.Stream):
        raise(Exception("Input transverse component is not a Stream object"))

    # Some integers
    nbin = len(RF_r)
    nt = len(RF_r[0].data)
    daz = 0.01
    naz = int(180./daz)
    deg2rad = np.pi/180.

    # Define time range over which to calculate azimuth
    indmin = int(t1/RF_r[0].stats.delta)
    indmax = int(t2/RF_r[0].stats.delta)

    # Copy stream stats
    str_stats = RF_r[0].stats

    # Initialize work arrays
    hr0_rot = np.zeros((nt, naz))
    ht0_rot = np.zeros((nt, naz))

    # Loop over each depth step
    for it in range(nt):

        # Initialize work arrays
        d_r = np.zeros(nbin)
        d_t = np.zeros(nbin)
        G = np.zeros((nbin, 5))

        # Build arrays and matrices 
        for itrace in range(nbin):

            baz = RF_r[itrace].stats.baz
            d_r[itrace] = RF_r[itrace].data[it]
            d_t[itrace] = RF_t[itrace].data[iz]
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
            hr0_rot[it, iaz] = np.cos(phi)*m_r[0]
            ht0_rot[it, iaz] = np.sin(phi)*m_t[0]


    # Minimize misfit of rotated transverse component over specific time range to
    # find azim
    RMS = np.zeros(naz)
    for iaz in range(naz):
        RMS[iaz] = np.sqrt(np.mean(np.square(ht0_rot[indmin:indmax, iaz])))

    # Azimuth of H1
    indaz = np.argmin(RMS)
    azcorr = indaz*daz

    # Resolve ambiguity based on radial component
    if hr0_rot[0,indaz]<0.:
        azcorr += 180.
        if azcorr > 360.:
            azcorr -= 360.

    return azcorr, RMS, hr0_rot, ht0_rot

