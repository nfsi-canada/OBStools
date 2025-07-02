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
import re
import math
import numpy as np
import fnmatch
from stdb import StDbElement
from scipy.stats import circmean
from matplotlib import pyplot as plt
from obspy import read, Stream, Trace, UTCDateTime
from obspy.core import AttribDict


def get_stkeys(inventory, keys=[]):

    allkeys = []
    station_list = inventory.get_contents()['stations']
    allkeys = [s.split(' ')[0] for s in station_list]

    if len(keys) > 0:
        # Extract key subset

        stkeys = []
        for key in keys:
            # Convert the pattern to a regex pattern
            # Replace '.' with '\.' to match literal dots
            # Replace '*' with '.*' to match any sequence of characters
            # Replace '?' with '.' to match any single character
            pattern = key.replace('.', r'\.').replace('*', '.*').replace('?', '.')

            # Compile the regex pattern
            regex = re.compile(f'^.*{pattern}.*$')

            # Filter allkeys based on the compiled regex
            stkeys.extend([key for key in allkeys if regex.match(key)])

    else:
        stkeys = allkeys

    return stkeys


def inv2stdb(inventory, keys=[]):

    stkeys = get_stkeys(inventory, keys)

    stations = {}
    for key in stkeys:
        net = key.split('.')[0]
        sta = key.split('.')[1]
        cha = '?H?'
        inv = inventory.select(network=net, station=sta, channel=cha)
        seed_id = inv.get_contents()['channels'][0]
        coords = inv.get_coordinates(seed_id)

        stdb_element = StDbElement(
            station=sta,
            network=net,
            channel=seed_id.split('.')[3][0:2],
            location=seed_id.split('.')[2],
            latitude=coords['latitude'],
            longitude=coords['longitude'],
            elevation=coords['elevation'],
            startdate=inv[0].stations[0].start_date,
            enddate=inv[0].stations[0].end_date
            )
        stations[key] = stdb_element

    return stations, stkeys


def traceshift(trace, tt):
    """
    Function to shift traces in time given travel time


    Parameters
    ----------

    trace : :class:`~obspy.core.Trace` object
        Trace object to update
    tt : float
        Time shift in seconds

    Returns
    -------

    rtrace : :class:`~obspy.core.Trace` object
        Updated trace object


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
    """
    Function for quality control of traces, which compares the
    start and end times that were requested, as well as the total n
    length of the traces.


    Parameters
    ----------

    start : :class:`~obspy.core.UTCDateTime` object
        Start time of requested stream
    end : :class:`~obspy.core.UTCDateTime` object
        End time of requested stream
    st : :class:`~obspy.core.Stream` object
        Stream object with all trace data

    Returns
    -------

    (pass): bool
        Whether the QC test has passed
    st : :class:`~obspy.core.Stream` object
        Updated stream object


    """

    # Check start times
    if not np.all([tr.stats.starttime == start for tr in st]):
        print("*      Start times are not all close to true start: ")
        [print("*        "+tr.stats.channel+" " +
               str(tr.stats.starttime)+" " +
               str(tr.stats.endtime)) for tr in st]
        print("*        True start: "+str(start))
        print("*   -> Shifting traces to true start")
        delay = [tr.stats.starttime - start for tr in st]
        st_shifted = Stream(
            traces=[traceshift(tr, dt) for tr, dt in zip(st, delay)])
        st = st_shifted.copy()

    # Try trimming
    dt = st[0].stats.delta
    try:
        st.trim(start, end-dt, fill_value=0., pad=True)
    except Exception:
        print("*   Unable to trim")
        print("*   -> Skipping")
        print("**************************************************")
        return False, None

    # Check final lengths - they should all be equal if start times
    # and sampling rates are all equal and traces have been trimmed
    sr = st[0].stats.sampling_rate
    if not np.allclose([tr.stats.npts for tr in st[1:]], st[0].stats.npts):
        print("*   Lengths are incompatible: ")
        [print("*       "+str(tr.stats.npts)) for tr in st]
        print("*   -> Skipping")
        print("**************************************************")

        return False, None

    elif not np.allclose([st[0].stats.npts], int((end - start)*sr), atol=1):
        print("*   Length is too short: ")
        print("*      "+str(st[0].stats.npts) +
              " ~= "+str(int((end - start)*sr)))
        print("*   -> Skipping")
        print("**************************************************")

        return False, None

    else:
        return True, st


def update_stats(tr, stla, stlo, stel, cha, evla=None, evlo=None):
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
    stel : float
        Station elevation (m)
    cha : str
        Channel for component
    evla : float, optional
        Latitude of event
    evlo : float, optional
        Longitute of event

    Returns
    -------

    tr : :class:`~obspy.core.Trace` object
        Updated trace object

    """

    tr.stats.sac = AttribDict()
    tr.stats.sac.stla = stla
    tr.stats.sac.stlo = stlo
    tr.stats.sac.stel = stel
    tr.stats.sac.kcmpnm = cha
    tr.stats.channel = cha
    if evla is not None and evlo is not None:
        tr.stats.sac.evla = evla
        tr.stats.sac.evlo = evlo

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
        Corresponding trace objects for components H1, H2, HZ and HP. Returns
        empty traces for missing components.

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
    Function to grab all available earthquake data given a path and data time
    range

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
        Corresponding trace objects for components H1, H2, HZ and HP. Returns
        empty traces for missing components.

    """

    # Find out how many events from Z.SAC files
    eventfiles = list(eventpath.glob('*Z.SAC'))
    if not eventfiles:
        raise Exception("No event found in folder "+str(eventpath))

    # Extract events from time stamps
    prefix = [file.name.split('.') for file in eventfiles]
    evstamp = [p[0]+'.'+p[1]+'.'+p[2]+'.'+p[3]+'.' for p in prefix]
    evDateTime = [UTCDateTime(p[0]+'-'+p[1]+'T'+p[2]+":"+p[3]) for p in prefix]

    # Define empty streams
    tr1 = Stream()
    tr2 = Stream()
    trZ = Stream()
    trP = Stream()

    # Cycle over all available files in time range
    for event, tstamp in zip(evDateTime, evstamp):
        if event >= tstart and event <= tend:

            # Cycle through directory and load files
            p = list(eventpath.glob('*.SAC'))
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


def calculate_tilt(ft1, ft2, ftZ, ftP, f, goodwins, tiltfreqs,
                   fig_trf=False, savefig=None):
    """
    Determines tilt orientation from the maximum coherence between
    rotated H1 and Z.

    Parameters
    ----------
    ft1, ft2, ftZ, ftP : :class:`~numpy.ndarray`
        Fourier transform of corresponding H1, H2, HZ and HP components
    f : :class:`~numpy.ndarray`
        Frequency axis in Hz
    goodwins : list
        List of booleans representing whether a window is good (True) or not
        (False). This attribute is returned from the method
        :func:`~obstools.atacr.classes.DayNoise.QC_daily_spectra`
    tiltfreqs : list
        Two floats representing the frequency band at which the mean
        coherence, phase and admittance are calculated to determine the
        tile orientation

    Returns
    -------
    cHH, cHZ, cHP : :class:`~numpy.ndarray`
        Arrays of power and cross-spectral density functions of components HH
        (rotated H1 in direction of maximum tilt), HZ, and HP
    coh : :class:`~numpy.ndarray`
        Mean coherence value between rotated H and Z components,
        as a function of azimuth CW from H1
    ph : :class:`~numpy.ndarray`
        Mean phase value between rotated H and Z components,
        as a function of azimuth CW from H1
    ad : :class:`~numpy.ndarray`
        Mean admittance value between rotated H and Z components,
        as a function of azimuth CW from H1
    phi : :class:`~numpy.ndarray`
        Array of azimuths CW from H1 considered
    tilt_dir : float
        Tilt direction (azimuth CW from H1) at maximum coherence
        between rotated H1 and Z
    tilt_ang : float
        Tilt angle (down from HZ/3) at maximum coherence
        between rotated H1 and Z
    coh_value : float
        Mean coherence value at tilt direction
    phase_value : float
        Mean phase value at tilt direction
    admit_value : float
        Mean admittance value at tilt direction

    """

    # frequencies considered in averaging
    freqs = (f > tiltfreqs[0]) & (f < tiltfreqs[1])

    # Mean of Z PSD - this one doesn't change with angle from H1
    cZZ = np.abs(np.mean(ftZ[goodwins, :] *
                         np.conj(ftZ[goodwins, :]), axis=0))[0:len(f)]

    # Initialize arrays
    phi = np.arange(0., 360., 10.)
    coh = np.zeros(len(phi))
    ph = np.zeros(len(phi))
    ad = np.zeros(len(phi))

    for i, d in enumerate(phi):

        # Rotate horizontals - clockwise from 1
        # components 1, 2 correspond to y, x in the cartesian plane
        ftH = rotate_dir(ft2, ft1, d)

        # Get transfer functions
        cHH = np.abs(np.mean(ftH[goodwins, :] *
                             np.conj(ftH[goodwins, :]), axis=0))[0:len(f)]
        cHZ = np.mean(ftH[goodwins, :] *
                      np.conj(ftZ[goodwins, :]), axis=0)[0:len(f)]

        Co = coherence(cHZ, cHH, cZZ)
        Ad = admittance(cHZ, cHH)
        Ph = phase(cHZ/cHH)

        # Calculate coherence over frequency band
        try:
            coh[i] = np.mean(Co[freqs])
            ad[i] = np.mean(Ad[freqs])
            ph[i] = circmean(
                Ph[freqs],
                high=np.pi,
                low=-np.pi)
        except Exception:
            print('Exception')
            coh[i] = 0.
            ad[i] = 0.
            ph[i] = 0.

    # Index where coherence is max
    ind = np.argwhere(coh == coh.max())

    # Phase, angle and direction at maximum coherence
    phase_value = ph[ind[0]][0]
    admit_value = ad[ind[0]][0]
    coh_value = coh[ind[0]][0]
    tilt_dir = phi[ind[0]][0]
    tilt_ang = np.arctan(admit_value)*180./np.pi

    # Phase has to be close to zero: ***Not true?***
    # Check phase near max coherence and 180 deg apart
    start = max(0, ind[0][0] - 1)
    end = min(len(ph), ind[0][0] + 2)
    phase_std = np.std(ph[start:end])

    if np.abs(phase_value) > np.pi/2:
        tilt_dir += 180.
    if tilt_dir > 360.:
        tilt_dir -= 360.

    # Refine search
    rphi = np.arange(tilt_dir-10., tilt_dir+10., 1.)
    rcoh = np.zeros(len(rphi))
    rad = np.zeros(len(rphi))
    rph = np.zeros(len(rphi))

    for i, d in enumerate(rphi):

        # Rotate horizontals - clockwise from 1
        # components 1, 2 correspond to y, x in the cartesian plane
        ftH = rotate_dir(ft2, ft1, d)

        # Get transfer functions
        cHH = np.abs(np.mean(ftH[goodwins, :] *
                             np.conj(ftH[goodwins, :]), axis=0))[0:len(f)]
        cHZ = np.mean(ftH[goodwins, :] *
                      np.conj(ftZ[goodwins, :]), axis=0)[0:len(f)]

        Co = coherence(cHZ, cHH, cZZ)
        Ad = admittance(cHZ, cHH)
        Ph = phase(cHZ/cHH)

        # Calculate coherence over frequency band
        try:
            rcoh[i] = np.mean(Co[freqs])
            rad[i] = np.mean(Ad[freqs])
            rph[i] = circmean(
                Ph[freqs],
                high=np.pi,
                low=-np.pi)
        except Exception:
            print("* Warning: problems in the Tilt calculations. " +
                  "Setting coherence and phase between Z and rotated H " +
                  "to 0.")
            rcoh[i] = 0.
            rad[i] = 0.
            rph[i] = 0.

    # Index where coherence is max
    ind = np.argwhere(rcoh == rcoh.max())

    # Phase and direction at maximum coherence
    phase_value = rph[ind[0]][0]
    admit_value = rad[ind[0]][0]
    coh_value = rcoh[ind[0]][0]
    tilt_dir = rphi[ind[0]][0]
    tilt_ang = np.arctan(admit_value)*180./np.pi

    # Now calculate spectra at tilt direction
    ftH = rotate_dir(ft2, ft1, tilt_dir)

    # Get transfer functions
    cHH = np.abs(np.mean(ftH[goodwins, :] *
                         np.conj(ftH[goodwins, :]), axis=0))[0:len(f)]
    cHZ = np.mean(ftH[goodwins, :] *
                  np.conj(ftZ[goodwins, :]), axis=0)[0:len(f)]
    if np.any(ftP):
        cHP = np.mean(ftH[goodwins, :] *
                      np.conj(ftP[goodwins, :]), axis=0)[0:len(f)]
    else:
        cHP = None

    return cHH, cHZ, cHP, coh, ph, ad, phi, tilt_dir, tilt_ang, coh_value, phase_value, admit_value


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
    axis : int, optional
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
    ----------
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
    ----------
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
    Calculates phase angle between real and imaginary components

    Parameters
    ----------
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


def rotate_dir(x, y, theta):
    """
    Rotates (x, y) data clockwise by an angle theta. Returns the
    rotated component y.

    Parameters
    ----------
    x : :class:`~numpy.ndarray`
        Array of values along the `x` coordinate
    y : :class:`~numpy.ndarray`
        Array of values along the `y` coordinate
    theta : :class:`~numpy.ndarray`
        Angle in degrees

    Returns
    -------
    y_rotated: :class:`~numpy.ndarray`
        Rotated array of component `y`

    """

    d = theta*np.pi/180.
    rot_mat = [[np.cos(d), np.sin(d)],
               [-np.sin(d), np.cos(d)]]

    vxy_rotated = np.tensordot(rot_mat, [x, y], axes=1)

    return vxy_rotated[1]


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


def robust(array):
    """
    Calculates robust quantities using the robust standard units.

    Parameters
    ----------
    array : :class:`~numpy.ndarray`
        Array of values

    Returns
    -------
    robust_array: :class:`~numpy.ndarray`
        Array with outliers removed
    outliers : :class:`~numpy.ndarray`
        Array of outlier values

    """

    median_array = np.median(array)
    mad_array = 1.4826*np.median(np.abs(array - median_array))
    if mad_array > 0.:
        rsu_array = (array - median_array)/mad_array
        robust_array = array[np.abs(rsu_array) < 2.]
        outliers = array[np.abs(rsu_array) >= 2.]
    else:
        robust_array = array
        outliers = None

    return robust_array, outliers
