'''
MODULE obs_util.py

Set of utility functions to calculate the transfer functions and perform
related operations

---
Pascal Audet
pascal.audet@uottawa.ca

Last updated: 17 November 2015

'''

import os
import numpy as np
import fnmatch
from matplotlib import pyplot as plt
from obspy.core import read, Stream, AttribDict
from scipy.signal import savgol_filter
from obs import obs_plot as obspl


def update_stats(tr, stla, stlo, stel):

    tr.stats.sac = AttribDict()
    tr.stats.sac.stla = stla
    tr.stats.sac.stlo = stlo
    tr.stats.sac.stel = stel
    
    return tr

def get_data(datapath, tstart, tend):

    # Define empty streams
    trN1 = Stream()
    trN2 = Stream()
    trNZ = Stream()
    trNP = Stream()

    # Time iterator
    t1 = tstart

    # Cycle through each day withing time range
    while t1 < tend:

        # Time stamp used in file name
        tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'

        # Cycle through directory and load files
        for file in os.listdir(datapath):
            if fnmatch.fnmatch(file, '*' + tstamp + '*1.SAC'):
                tr = read(datapath + file)
                trN1.append(tr[0])
            elif fnmatch.fnmatch(file, '*' + tstamp + '*2.SAC'):
                tr = read(datapath + file)
                trN2.append(tr[0])
            elif fnmatch.fnmatch(file, '*' + tstamp + '*Z.SAC'):
                tr = read(datapath + file)
                trNZ.append(tr[0])
            elif fnmatch.fnmatch(file, '*' + tstamp + '*H.SAC'):
                tr = read(datapath + file)
                trNP.append(tr[0])

        # Increase increment
        t1 += 3600.*24.

    return trN1, trN2, trNZ, trNP


def get_event(eventpath, tstart, tend):

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
        for file in os.listdir(datapath):
            if fnmatch.fnmatch(file, '*' + tstamp + '*1.SAC'):
                tr = read(datapath + file)
                tr1.append(tr[0])
            elif fnmatch.fnmatch(file, '*' + tstamp + '*2.SAC'):
                tr = read(datapath + file)
                tr2.append(tr[0])
            elif fnmatch.fnmatch(file, '*' + tstamp + '*Z.SAC'):
                tr = read(datapath + file)
                trZ.append(tr[0])
            elif fnmatch.fnmatch(file, '*' + tstamp + '*H.SAC'):
                tr = read(datapath + file)
                trP.append(tr[0])

        # Increase increment
        t1 += 3600.*24.

    return trN1, trN2, trNZ, trNP

def calculate_tilt(ft1, ft2, ftZ, ftP, f, goodwins, tiltfreq=[0.005, 0.035]):
    """ 
    Determines tilt direction from maximum
    coherence between H1 and Z

    """

    direc = np.arange(0., 360., 10.)
    coh = np.zeros(len(direc))
    ph = np.zeros(len(direc))
    cZZ = np.abs(np.mean(ftZ[goodwins,:]*np.conj(ftZ[goodwins,:]), axis=0))[0:len(f)]

    for i, d in enumerate(direc):

        # Rotate horizontals
        ftH = rotate_dir(ft1, ft2, d)
    
        # Get transfer functions
        cHH = np.abs(np.mean(ftH[goodwins,:]*np.conj(ftH[goodwins,:]), axis=0))[0:len(f)]
        cHZ = np.mean(ftH[goodwins,:]*np.conj(ftZ[goodwins,:]), axis=0)[0:len(f)]

        Co = coherence(cHZ, cHH, cZZ)
        Ph = phase(cHZ)

        # Calculate coherence over frequency band
        coh[i] = np.mean(Co[(f>tiltfreq[0]) & (f<tiltfreq[1])])
        ph[i] = np.pi/2. - np.mean(Ph[(f>tiltfreq[0]) & (f<tiltfreq[1])])

    # Index where coherence is max
    ind = np.argwhere(coh==coh.max())

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
        cHH = np.abs(np.mean(ftH[goodwins,:]*np.conj(ftH[goodwins,:]), axis=0))[0:len(f)]
        cHZ = np.mean(ftH[goodwins,:]*np.conj(ftZ[goodwins,:]), axis=0)[0:len(f)]

        Co = coherence(cHZ, cHH, cZZ)
        Ph = phase(cHZ)

        # Calculate coherence over frequency band
        rcoh[i] = np.mean(Co[(f>tiltfreq[0]) & (f<tiltfreq[1])])
        rph[i] = np.pi/2. - np.mean(Ph[(f>tiltfreq[0]) & (f<tiltfreq[1])])

    # Index where coherence is max
    ind = np.argwhere(rcoh==rcoh.max())

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
    cHH = np.abs(np.mean(ftH[goodwins,:]*np.conj(ftH[goodwins,:]), axis=0))[0:len(f)]
    cHZ = np.mean(ftH[goodwins,:]*np.conj(ftZ[goodwins,:]), axis=0)[0:len(f)]
    cHP = np.mean(ftH[goodwins,:]*np.conj(ftP[goodwins,:]), axis=0)[0:len(f)]

    return cHH, cHZ, cHP, coh, ph, direc, tilt, coh_value, phase_value


def calculate_windowed_fft(trace, ws, ss=None, hann=True):
    """
    Calculates cross-spectral quantities

    ws: window size
    ss: step size (number of samples until next window)

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
    if data.ndim > 1:
        filt = np.zeros(data.shape)
        for i in range(data.shape[::-1][axis]):
            if axis==0:
                filt[:,i] = np.convolve(data[:,i], np.ones((nd,))/nd, mode='same')
            elif axis==1:
                filt[i,:] = np.convolve(data[i,:], np.ones((nd,))/nd, mode='same')
    else:
        filt = np.convolve(data, np.ones((nd,))/nd, mode='same')
    return filt


def admittance(Gxy, Gxx):

    return np.abs(Gxy)/Gxx


def coherence(Gxy, Gxx, Gyy):

    return np.abs(Gxy)**2/(Gxx*Gyy)


def phase(Gxy):
    
    return np.angle(Gxy)


def sliding_window(a, ws, ss=None, hann=True):
    '''
    Parameters
        a  - a 1D array
        ws - the window size, in samples
        ss - the step size, in samples. If not provided, window and step size
             are equal.

    """ PA: This is not my function """

    '''
     
    if ss is None:
        # no step size was provided. Return non-overlapping windows
        ss = ws
    
    # Calculate the number of windows to return, ignoring leftover samples, and
    # allocate memory to contain the samples
    valid = len(a) - ss
    nd = (valid) // ss
    out = np.ndarray((nd,ws), dtype=a.dtype)

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
            out[i] = a[start : stop] * np.hanning(ws)
        else:
            out[i] = a[start : stop]
     
    return out, nd


def rotate_dir(tr1, tr2, direc):
    """
    Rotate horizontals

    """

    d = -direc*np.pi/180.+np.pi/2.
    rot_mat = np.array([[np.cos(d), -np.sin(d)],
        [np.sin(d), np.cos(d)]])

    v12 = np.array([tr2, tr1])
    vxy = np.tensordot(rot_mat, v12, axes=1)
    tr_2 = vxy[0,:]
    tr_1 = vxy[1,:]

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
    return 1 if x==0 else 2**(x-1).bit_length()

# def get_trZ_p(tr, pAd, pPh, f1, f2):
#     """ 
#     Calculates predicted trace from transfer 
#     functions

#     """
#     # Copy trace to predicted vertical (for consistency)
#     tr_p = tr.copy()
    
#     # Fourier transform trace 
#     ftr = np.fft.fft(tr.data)

#     nt = tr.stats.npts
#     dt = tr.stats.delta
#     ff = np.fft.fftfreq(nt, dt)

#     Ad = np.zeros((nt,))
#     Ph = np.zeros((nt,))
#     ind1 = np.where((ff>f1)&(ff<f2))
#     ind2 = np.where((ff<-f1)&(ff>-f2))

#     # Define analytical transfer functions
#     # Careful: phase of negative frequencies is negative
#     Ad[ind1] = pAd[0]*ff[ind1]**2 + pAd[1]*ff[ind1] + pAd[2]
#     Ad[ind2] = pAd[0]*ff[ind2]**2 + pAd[1]*ff[ind2] + pAd[2]
#     Ph[ind1] = pPh[0]*ff[ind1]**2 + pPh[1]*ff[ind1] + pPh[2]
#     Ph[ind2] = -(pPh[0]*ff[ind2]**2 + pPh[1]*ff[ind2] + pPh[2])

#     # Define complex transfer function in frequency
#     trf = Ad*np.exp(1j*Ph)

#     # Convolve in Fourier domain
#     ftr_p = trf*ftr

#     # Inverse Fourier transform
#     tr_p.data = np.real(np.fft.ifft(ftr_p))

#     return tr_p


# def calculate_xspecs(data1, data2, goodwins):
#     """
#     Calculates cross-spectral quantities

#     ws: window size
#     ss: step size (number of samples until next window)

#     """

#     # data1 = data1[]
#     Gxx = np.abs(np.mean(np.conj(ftx)*ftx, axis=0))
#     Gyy = np.abs(np.mean(np.conj(fty)*fty, axis=0))
#     Gxy = np.mean(np.conj(ftx)*fty, axis=0)

#     Cxy = np.mean(np.real(ftx)*np.real(fty) + 
#             np.imag(ftx)*np.imag(fty), axis=0)
#     Qxy = np.mean(np.real(ftx)*np.imag(fty) - 
#             np.imag(ftx)*np.real(fty), axis=0)

#     return Gxy, Gxx, Gyy, Cxy, Qxy, nd


# def get_trfs_lsqfit(Ad, Co, Ph, eAd, eCo, ePh, freq, f1, f2, plot=True):
#     """
#     Get best-fit analytical admittance and phase from observed
#     admittance and phase data within a given frequency range
    
#     """

#     ind = np.where((freq>f1)&(freq<f2))
#     cond = np.where((Co>0.5)&(freq>f1)&(freq<f2))
#     if not cond:
#         print('no coherence above 0.5 within ',f1, f2, 'Hz')
#         return [0., 0., 0.], [0., 0., 0.]    
#     mcoh = np.mean(Co[ind])
#     if mcoh < 0.5:
#         print('mean coherence too low', mcoh)
#         return [0., 0., 0.], [0., 0., 0.]


#     pAd = np.polyfit(freq[cond],Ad[cond],2,w=1./eAd[cond])
#     pPh = np.polyfit(freq[cond],Ph[cond],2,w=1./ePh[cond])

#     aAd = pAd[0]*freq**2 + pAd[1]*freq + pAd[2]
#     aPh = pPh[0]*freq**2 + pPh[1]*freq + pPh[2]

#     if plot:
#         obspl.obs_plot_trf_analytic(Ad, Ph, eAd, ePh, aAd, aPh, freq, ind)

#     return pAd, pPh



# def trfs_error(Ad, Co, Ph, nd):
#     """
#     Errors on transfer functions

#     """

#     eps = np.sqrt((1.-Co)/(2.*Co*nd))

#     eAd = Ad*eps
#     eCo = Co*eps
#     ePh = eps

#     return eAd, eCo, ePh


