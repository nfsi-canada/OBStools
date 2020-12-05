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
:mod:`~obstools.atacr.utils` contains several functions that are used in the
class methods of `~obstools.atacr.classes`.

"""


import os
import math
import numpy as np
import fnmatch
from matplotlib import pyplot as plt
from obspy.core import read, Stream, Trace, AttribDict
from scipy.signal import savgol_filter


def floor_decimal(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier) / multiplier


def traceshift(trace, tt):
    """
    Function to shift traces in time given travel time

    """

    # Define frequencies
    nt = trace.stats.npts
    dt = trace.stats.delta
    freq = np.fft.fftfreq(nt, d=dt)

    # Fourier transform
    ftrace = np.fft.fft(trace.data)

    # Shift
    for i in range(len(freq)):
        ftrace[i] = ftrace[i]*np.exp(-2.*np.pi*1j*freq[i]*tt)

    # Back Fourier transform and return as trace
    rtrace = trace.copy()
    rtrace.data = np.real(np.fft.ifft(ftrace))

    # Update start time
    rtrace.stats.starttime -= tt

    return rtrace


def QC_streams(start, end, st):

    # Check start times
    if not np.all([tr.stats.starttime == start for tr in st]):
        print("* Start times are not all close to true start: ")
        [print("*   "+tr.stats.channel+" " +
               str(tr.stats.starttime)+" " +
               str(tr.stats.endtime)) for tr in st]
        print("*   True start: "+str(start))
        print("* -> Shifting traces to true start")
        delay = [tr.stats.starttime - start for tr in st]
        st_shifted = Stream(
            traces=[traceshift(tr, dt) for tr, dt in zip(st, delay)])
        st = st_shifted.copy()

    # # Check sampling rate
    # sr = st[0].stats.sampling_rate
    # sr_round = float(floor_decimal(sr, 0))
    # if not sr == sr_round:
    #     print("* Sampling rate is not an integer value: ", sr)
    #     print("* -> Resampling")
    #     st.resample(sr_round, no_filter=False)

    # Try trimming
    dt = st[0].stats.delta
    try:
        st.trim(start, end-dt, fill_value=0., pad=True)
    except:
        print("* Unable to trim")
        print("* -> Skipping")
        print("**************************************************")
        return False, None

    # Check final lengths - they should all be equal if start times
    # and sampling rates are all equal and traces have been trimmed
    sr = st[0].stats.sampling_rate
    if not np.allclose([tr.stats.npts for tr in st[1:]], st[0].stats.npts):
        print("* Lengths are incompatible: ")
        [print("*     "+str(tr.stats.npts)) for tr in st]
        print("* -> Skipping")
        print("**************************************************")

        return False, None

    elif not np.allclose([st[0].stats.npts], int((end - start)*sr), atol=1):
        print("* Length is too short: ")
        print("*    "+str(st[0].stats.npts) +
              " ~= "+str(int((end - start)*sr)))
        print("* -> Skipping")
        print("**************************************************")

        return False, None

    else:
        return True, st


def update_stats(tr, stla, stlo, stel, cha):
    """
    Function to include SAC metadata to :class:`~obspy.core.Trace` objects

    Parameters
    ----------

    tr : :class:`~obspy.core.Trace` object
        Trace object to update
    stla : float
        Latitude of station
    stlo : float
        Longitude of station
    cha : str
        Channel for component

    Returns
    -------

    tr : :class:`~obspy.core.Trace` object
        Updated trace object

    """

    tr.stats.sac = AttribDict()
    tr.stats.sac.stla = stla
    tr.stats.sac.stlo = stlo
    tr.stats.sac.stel = stel
    tr.stats.channel = cha

    return tr


def get_data(datapath, tstart, tend):
    """
    Function to grab all available noise data given a path and data time range

    Parameters
    ----------
    datapath : str
        Path to noise data folder
    tstart : :class:`~obspy.class.UTCDateTime`
        Start time for query
    tend : :class:`~obspy.class.UTCDateTime`
        End time for query

    Returns
    -------
    tr1, tr2, trZ, trP : :class:`~obspy.core.Trace` object
        Corresponding trace objects for components H1, H2, HZ and HP. Returns empty traces
        for missing components.

    """

    # Define empty streams
    trN1 = Stream()
    trN2 = Stream()
    trNZ = Stream()
    trNP = Stream()

    # Time iterator
    t1 = tstart

    # Cycle through each day within time range
    while t1 < tend:

        # Time stamp used in file name
        tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'

        # Cycle through directory and load files
        p = datapath.glob('*.*')
        files = [x for x in p if x.is_file()]
        for file in files:
            if fnmatch.fnmatch(str(file), '*' + tstamp + '*1.SAC'):
                tr = read(str(file))
                trN1.append(tr[0])
            elif fnmatch.fnmatch(str(file), '*' + tstamp + '*2.SAC'):
                tr = read(str(file))
                trN2.append(tr[0])
            elif fnmatch.fnmatch(str(file), '*' + tstamp + '*Z.SAC'):
                tr = read(str(file))
                trNZ.append(tr[0])
            elif fnmatch.fnmatch(str(file), '*' + tstamp + '*H.SAC'):
                tr = read(str(file))
                trNP.append(tr[0])

        # Increase increment
        t1 += 3600.*24.

    # Fill with empty traces if components are not found
    ntr = len(trNZ)
    if not trN1 and not trN2:
        for i in range(ntr):
            trN1.append(Trace())
            trN2.append(Trace())
    if not trNP:
        for i in range(ntr):
            trNP.append(Trace())

    if ntr > 0:
        # Check that all sampling rates are equal - otherwise resample
        if trNZ[0].stats.sampling_rate != trNP[0].stats.sampling_rate:

            # These checks assume that all seismic data have the same sampling
            if trNZ[0].stats.sampling_rate < trNP[0].stats.sampling_rate:
                trNP.resample(trNZ[0].stats.sampling_rate, no_filter=False)
            else:
                trNZ.resample(trNP[0].stats.sampling_rate, no_filter=False)
                if trN1:
                    trN1.resample(trNP[0].stats.sampling_rate, no_filter=False)
                if trN2:
                    trN2.resample(trNP[0].stats.sampling_rate, no_filter=False)

    return trN1, trN2, trNZ, trNP


def get_event(eventpath, tstart, tend):
    """
    Function to grab all available earthquake data given a path and data time range

    Parameters
    ----------
    eventpath : str
        Path to earthquake data folder
    tstart : :class:`~obspy.class.UTCDateTime`
        Start time for query
    tend : :class:`~obspy.class.UTCDateTime`
        End time for query

    Returns
    -------
    tr1, tr2, trZ, trP : :class:`~obspy.core.Trace` object
        Corresponding trace objects for components H1, H2, HZ and HP. Returns empty traces
        for missing components.

    """

    # Define empty streams
    tr1 = Stream()
    tr2 = Stream()
    trZ = Stream()
    trP = Stream()

    # Time iterator
    t1 = tstart

    # Cycle through each day withing time range
    while t1 < tend:

        # Time stamp used in file name
        tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'

        # Cycle through directory and load files
        p = eventpath.glob('*.*')
        files = [x for x in p if x.is_file()]
        for file in files:
            if fnmatch.fnmatch(str(file), '*' + tstamp + '*1.SAC'):
                tr = read(str(file))
                tr1.append(tr[0])
            elif fnmatch.fnmatch(str(file), '*' + tstamp + '*2.SAC'):
                tr = read(str(file))
                tr2.append(tr[0])
            elif fnmatch.fnmatch(str(file), '*' + tstamp + '*Z.SAC'):
                tr = read(str(file))
                trZ.append(tr[0])
            elif fnmatch.fnmatch(str(file), '*' + tstamp + '*H.SAC'):
                tr = read(str(file))
                trP.append(tr[0])

        # Increase increment
        t1 += 3600.*24.

    # Fill with empty traces if components are not found
    ntr = len(trZ)
    if not tr1 and not tr2:
        for i in range(ntr):
            tr1.append(Trace())
            tr2.append(Trace())
    if not trP:
        for i in range(ntr):
            trP.append(Trace())

    if ntr > 0:
        # Check that all sampling rates are equal - otherwise resample
        if trZ[0].stats.sampling_rate != trP[0].stats.sampling_rate:

            # These checks assume that all seismic data have the same sampling
            if trZ[0].stats.sampling_rate < trP[0].stats.sampling_rate:
                trP.resample(trZ[0].stats.sampling_rate, no_filter=False)
            else:
                trZ.resample(trP[0].stats.sampling_rate, no_filter=False)
                if tr1:
                    tr1.resample(trP[0].stats.sampling_rate, no_filter=False)
                if tr2:
                    tr2.resample(trP[0].stats.sampling_rate, no_filter=False)

    return tr1, tr2, trZ, trP


def calculate_tilt(ft1, ft2, ftZ, ftP, f, goodwins, tiltfreq=[0.005, 0.035]):
    """
    Determines tilt direction from maximum coherence between rotated H1 and Z.

    Parameters
    ----------
    ft1, ft2, ftZ, ftP : :class:`~numpy.ndarray`
        Fourier transform of corresponding H1, H2, HZ and HP components
    f : :class:`~numpy.ndarray`
        Frequency axis in Hz
    goodwins : list
        List of booleans representing whether a window is good (True) or not (False).
        This attribute is returned from the method :func:`~obstools.atacr.classes.DayNoise.QC_daily_spectra`
    tiltfreq : list
        Two floats representing the frequency band at which the tilt is calculated

    Returns
    -------
    cHH, cHZ, cHP : :class:`~numpy.ndarray`
        Arrays of power and cross-spectral density functions of components HH (rotated H1
        in direction of maximum tilt), HZ, and HP
    coh : :class:`~numpy.ndarray`
        Coherence value between rotated H and Z components, as a function of directions (azimuths)
    ph : :class:`~numpy.ndarray`
        Phase value between rotated H and Z components, as a function of directions (azimuths)
    direc : :class:`~numpy.ndarray`
        Array of directions (azimuths) considered
    tilt : float
        Direction (azimuth) of maximum coherence between rotated H1 and Z
    coh_value : float
        Coherence value at tilt direction
    phase_value : float
        Phase value at tilt direction

    """

    direc = np.arange(0., 360., 10.)
    coh = np.zeros(len(direc))
    ph = np.zeros(len(direc))
    cZZ = np.abs(np.mean(ftZ[goodwins, :] *
                         np.conj(ftZ[goodwins, :]), axis=0))[0:len(f)]

    for i, d in enumerate(direc):

        # Rotate horizontals
        ftH = rotate_dir(ft1, ft2, d)

        # Get transfer functions
        cHH = np.abs(np.mean(ftH[goodwins, :] *
                             np.conj(ftH[goodwins, :]), axis=0))[0:len(f)]
        cHZ = np.mean(ftH[goodwins, :] *
                      np.conj(ftZ[goodwins, :]), axis=0)[0:len(f)]

        Co = coherence(cHZ, cHH, cZZ)
        Ph = phase(cHZ)

        # Calculate coherence over frequency band
        coh[i] = np.mean(Co[(f > tiltfreq[0]) & (f < tiltfreq[1])])
        ph[i] = np.pi/2. - np.mean(Ph[(f > tiltfreq[0]) & (f < tiltfreq[1])])

    # Index where coherence is max
    ind = np.argwhere(coh == coh.max())

    # Phase and direction at maximum coherence
    phase_value = ph[ind[0]][0]
    coh_value = coh[ind[0]][0]
    tilt = direc[ind[0]][0]

    # Refine search
    rdirec = np.arange(direc[ind[0]][0]-10., direc[ind[0]][0]+10., 1.)
    rcoh = np.zeros(len(direc))
    rph = np.zeros(len(direc))

    for i, d in enumerate(rdirec):

        # Rotate horizontals
        ftH = rotate_dir(ft1, ft2, d)

        # Get transfer functions
        cHH = np.abs(np.mean(ftH[goodwins, :] *
                             np.conj(ftH[goodwins, :]), axis=0))[0:len(f)]
        cHZ = np.mean(ftH[goodwins, :] *
                      np.conj(ftZ[goodwins, :]), axis=0)[0:len(f)]

        Co = coherence(cHZ, cHH, cZZ)
        Ph = phase(cHZ)

        # Calculate coherence over frequency band
        rcoh[i] = np.mean(Co[(f > tiltfreq[0]) & (f < tiltfreq[1])])
        rph[i] = np.pi/2. - np.mean(Ph[(f > tiltfreq[0]) & (f < tiltfreq[1])])

    # Index where coherence is max
    ind = np.argwhere(rcoh == rcoh.max())

    # Phase and direction at maximum coherence
    phase_value = rph[ind[0]][0]
    coh_value = rcoh[ind[0]][0]
    tilt = rdirec[ind[0]][0]

    # Phase has to be close to zero - otherwise add pi
    if phase_value > 0.5*np.pi:
        tilt += 180.
    if tilt > 360.:
        tilt -= 360.

    # print('Maximum coherence for tilt = ', tilt)

    # Now calculate spectra at tilt direction
    ftH = rotate_dir(ft1, ft2, tilt)

    # Get transfer functions
    cHH = np.abs(np.mean(ftH[goodwins, :] *
                         np.conj(ftH[goodwins, :]), axis=0))[0:len(f)]
    cHZ = np.mean(ftH[goodwins, :]*np.conj(ftZ[goodwins, :]), axis=0)[0:len(f)]
    if np.any(ftP):
        cHP = np.mean(ftH[goodwins, :] *
                      np.conj(ftP[goodwins, :]), axis=0)[0:len(f)]
    else:
        cHP = None

    return cHH, cHZ, cHP, coh, ph, direc, tilt, coh_value, phase_value


def calculate_windowed_fft(trace, ws, ss=None, hann=True):
    """
    Calculates windowed Fourier transform

    Parameters
    ----------
    trace : :class:`~obspy.core.Trace`
        Input trace data
    ws : int
        Window size, in number of samples
    ss : int
        Step size, or number of samples until next window
    han : bool
        Whether or not to apply a Hanning taper to data

    Returns
    -------
    ft : :class:`~numpy.ndarray`
        Fourier transform of trace
    f : :class:`~numpy.ndarray`
        Frequency axis in Hz

    """

    n2 = _npow2(ws)
    f = trace.stats.sampling_rate/2. * np.linspace(0., 1., int(n2/2) + 1)

    # Extract sliding windows
    tr, nd = sliding_window(trace.data, ws, ss)

    # Fourier transform
    ft = np.fft.fft(tr, n=n2)

    return ft, f

# def smooth(data, np, poly=0, axis=0):
#     return savgol_filter(data, np, poly, axis=axis, mode='wrap')


def smooth(data, nd, axis=0):
    """
    Function to smooth power spectral density functions from the convolution
    of a boxcar function with the PSD

    Parameters
    ----------
    data : :class:`~numpy.ndarray`
        Real-valued array to smooth (PSD)
    nd : int
        Number of samples over which to smooth
    axis : int
        axis over which to perform the smoothing

    Returns
    -------
    filt : :class:`~numpy.ndarray`, optional
        Filtered data

    """
    if np.any(data):
        if data.ndim > 1:
            filt = np.zeros(data.shape)
            for i in range(data.shape[::-1][axis]):
                if axis == 0:
                    filt[:, i] = np.convolve(
                        data[:, i], np.ones((nd,))/nd, mode='same')
                elif axis == 1:
                    filt[i, :] = np.convolve(
                        data[i, :], np.ones((nd,))/nd, mode='same')
        else:
            filt = np.convolve(data, np.ones((nd,))/nd, mode='same')
        return filt
    else:
        return None


def admittance(Gxy, Gxx):
    """
    Calculates admittance between two components

    Parameters
    ---------
    Gxy : :class:`~numpy.ndarray`
        Cross spectral density function of `x` and `y`
    Gxx : :class:`~numpy.ndarray`
        Power spectral density function of `x`

    Returns
    -------
    : :class:`~numpy.ndarray`, optional
        Admittance between `x` and `y`

    """

    if np.any(Gxy) and np.any(Gxx):
        return np.abs(Gxy)/Gxx
    else:
        return None


def coherence(Gxy, Gxx, Gyy):
    """
    Calculates coherence between two components

    Parameters
    ---------
    Gxy : :class:`~numpy.ndarray`
        Cross spectral density function of `x` and `y`
    Gxx : :class:`~numpy.ndarray`
        Power spectral density function of `x`
    Gyy : :class:`~numpy.ndarray`
        Power spectral density function of `y`

    Returns
    -------
    : :class:`~numpy.ndarray`, optional
        Coherence between `x` and `y`

    """

    if np.any(Gxy) and np.any(Gxx) and np.any(Gxx):
        return np.abs(Gxy)**2/(Gxx*Gyy)
    else:
        return None


def phase(Gxy):
    """
    Calculates phase angle between two components

    Parameters
    ---------
    Gxy : :class:`~numpy.ndarray`
        Cross spectral density function of `x` and `y`

    Returns
    -------
    : :class:`~numpy.ndarray`, optional
        Phase angle between `x` and `y`

    """

    if np.any(Gxy):
        return np.angle(Gxy)
    else:
        return None


def sliding_window(a, ws, ss=None, hann=True):
    """
    Function to split a data array into overlapping, possibly tapered sub-windows

    Parameters
    ----------
    a : :class:`~numpy.ndarray`
        1D array of data to split
    ws : int
        Window size in samples
    ss : int
        Step size in samples. If not provided, window and step size
         are equal.

    Returns
    -------
    out : :class:`~numpy.ndarray`
        1D array of windowed data
    nd : int
        Number of windows

    """

    if ss is None:
        # no step size was provided. Return non-overlapping windows
        ss = ws

    # Calculate the number of windows to return, ignoring leftover samples, and
    # allocate memory to contain the samples
    valid = len(a) - ss
    nd = (valid) // ss
    out = np.ndarray((nd, ws), dtype=a.dtype)

    if nd == 0:
        if hann:
            out = a * np.hanning(ws)
        else:
            out = a

    for i in range(nd):
        # "slide" the window along the samples
        start = i * ss
        stop = start + ws
        if hann:
            out[i] = a[start: stop] * np.hanning(ws)
        else:
            out[i] = a[start: stop]

    return out, nd


def rotate_dir(tr1, tr2, direc):

    d = -direc*np.pi/180.+np.pi/2.
    rot_mat = np.array([[np.cos(d), -np.sin(d)],
                        [np.sin(d), np.cos(d)]])

    v12 = np.array([tr2, tr1])
    vxy = np.tensordot(rot_mat, v12, axes=1)
    tr_2 = vxy[0, :]
    tr_1 = vxy[1, :]

    return tr_1


def ftest(res1, pars1, res2, pars2):

    from scipy.stats import f as f_dist

    N1 = len(res1)
    N2 = len(res2)

    dof1 = N1 - pars1
    dof2 = N2 - pars2

    Ea_1 = np.sum(res1**2)
    Ea_2 = np.sum(res2**2)

    Fobs = (Ea_1/dof1)/(Ea_2/dof2)

    P = 1. - (f_dist.cdf(Fobs, dof1, dof2) - f_dist.cdf(1./Fobs, dof1, dof2))

    return P


def _npow2(x):
    return 1 if x == 0 else 2**(x-1).bit_length()
