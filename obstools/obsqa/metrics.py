from .metrics import *
import numpy as _np
# from scipy.signal import csd as _csd
import obspy
from scipy.signal import stft, detrend

def _window(dt,window=None,overlap=0.3):
        # Points in window
        ws = int(window/dt)
        # Number of points to overlap
        ss = int(window*overlap/dt)
        # hanning window
        hanning = _np.hanning(2*ss)
        wind = _np.ones(ws)
        wind[0:ss] = hanning[0:ss]
        wind[-ss:ws] = hanning[ss:ws]
        return wind,ws,ss
def _stft(tr,fs,overlap=0.3,window=None,return_onesided=False,scaling='spectrum'):
        if window==None:
                tlen = _np.size(tr)/fs>7200
                if tlen:
                        window = 7200
                else:
                        window = 500
        wind,ws,ss = _window(1/fs,window=window,overlap=overlap)
        _f, _t, ft = stft(tr,fs, return_onesided=return_onesided, boundary=None,padded=False, window=wind, nperseg=ws, noverlap=ss,detrend='constant',scaling=scaling)
        ft = ft*ws
        ft = _np.mean(ft,axis=1)
        return _f.T,ft.T
def _psd(tr,fs,window=None,overlap=0.3):
        f,ft = _stft(tr,fs,overlap=overlap,scaling='psd',return_onesided=True,window=window)
        psd = _np.abs(ft)**2
        return f,psd
def _csd(a_ft,b_ft):
        cab = a_ft*_np.conj(b_ft)
        return cab
def _csd_helper(a,b,fs,return_onesided=False,window=None):
        f,a_ft = _stft(a,fs,return_onesided=return_onesided,window=window)
        f,b_ft = _stft(b,fs,return_onesided=return_onesided,window=window)
        cab = _csd(a_ft,b_ft)
        return f,cab

def Phase(ab):
        '''
        Takes a single NDarray containing the cross-psd (complex) between two components and pulls the phase out of the imaginary component of the array. Output is in degrees for a quadrature circle (-180,180).
        '''
        ph = _np.angle(ab,deg=True)
        return ph

def Admittance(ab,bb):
        '''
        Takes two NDarrays of equal shape, the first (a) containing the cross-power-spectral density between the two components and the second (b) containing the auto-power-spectral density of the secondary component.
        ie, For the spectral admittance between Z and P, a is the cross psd between the two and b is the auto-psd of P.
        '''
        ad = _np.abs(ab)/bb
        return ad
def Coherence(ab,aa,bb):
        '''
        Takes three NDarrays of equal shape, the first (ab) containing the cross-power-spectral density between the two components, the second and third (aa and bb) containing the auto-power-spectral density of these two components.
        ie, For the spectral coherence between Z and P, a is the cross psd between the two and b and c are the auto-psd of Z and P, respectively.
        '''
        coh = ((abs(ab)**2)/(aa*bb))
        return coh
def _csd(a_ft,b_ft):
        cab = a_ft*_np.conj(b_ft)
        return cab
def _csd_helper(a,b,fs,return_onesided=False,window=None):
        f,a_ft = _stft(a,fs,return_onesided=return_onesided,window=window)
        f,b_ft = _stft(b,fs,return_onesided=return_onesided,window=window)
        cab = _csd(a_ft,b_ft)
        return f,cab

def CrossSpec(A,B,fs=None,return_onesided=False,window=None):
        if isinstance(A,obspy.core.trace.Trace):
                fs = 1/A.stats.delta
                A = A.data
        if isinstance(B,obspy.core.trace.Trace):
                fs = 1/B.stats.delta
                B = B.data
        f,spec_AB= _csd_helper(A,B, fs=fs,return_onesided=return_onesided,window=window)
        f,spec_AA= _csd_helper(A,A, fs=fs,return_onesided=return_onesided,window=window)
        f,spec_BB= _csd_helper(B,B, fs=fs,return_onesided=return_onesided,window=window)
        return f,spec_AB,_np.abs(spec_AA),_np.abs(spec_BB)

def CrossCOH(A,B,fs=None):
        f,spec_AB,spec_AA,spec_BB = CrossSpec(A,B,fs=fs)
        coh = Coherence(spec_AB,spec_AA,spec_BB)
        coh = coh[f>=0]
        f = f[f>=0]
        return f,coh

def CrossPh(A,B,fs=None):
        f,spec_AB,spec_AA,spec_BB = CrossSpec(A,B,fs=fs)
        ph = Phase(spec_AB)
        return f,ph

def CrossAdm(A,B,fs=None):
        f,spec_AB,spec_AA,spec_BB = CrossSpec(A,B,fs=fs)
        adm = Admittance(spec_AB,spec_BB)
        return f,adm

def Metrics(A,B,fs=None,return_onesided=False,window=None):
        f,spec_AB,spec_AA,spec_BB = CrossSpec(A,B,fs=fs,return_onesided=return_onesided,window=window)
        coh = Coherence(spec_AB,spec_AA,spec_BB)
        adm = Admittance(spec_AB,spec_BB)
        ph = Phase(spec_AB)
        return f,coh,adm,ph

def obmetrics(spectra,meters=['Coherence','Admittance','Phase'],crosskeys = ['12','1Z','2Z','ZP']):
        # Takes a given dictionary containing cross (& auto) spectral densities
        # and returns their cross coherence, phase, and admittance.
        try:
                f = spectra['f']
        except:
                f = []
        result = dict()
        for mtr in meters:
                result[mtr] = dict()
                for ckey in crosskeys:
                        ab = spectra['c' + ckey]
                        aa = spectra['c' + ckey[0] + ckey[0]]
                        bb = spectra['c' + ckey[1] + ckey[1]]
                        if mtr=='Coherence':
                                result[mtr][ckey] = Coherence(ab,aa,bb)
                        if mtr=='Phase':
                                result[mtr][ckey] = Phase(ab)
                        if mtr=='Admittance':
                                result[mtr][ckey] = Admittance(ab,bb)
                        if len(f)>0:
                                result[mtr][ckey] = result[mtr][ckey].reshape(-1)[f.reshape(-1)>=0]
        if len(f)>0:
                f = f.reshape(-1)
                f = f[f>=0]
        result['f'] = f
        return result

def rad(x,N=20,logfilt=False):
    '''Rolling Absolute Deviation (Spike Detector) - Charles H. '23
    logfilt: Filters to zero all values not more than a magnitude larger than the average.
    This usually removes the noisefloor entirely but will kill any signal with an SNR<10.'''
    out = _np.abs(x - _np.abs(_np.convolve(x, _np.ones(int((N)))/int((N)), mode='same')))
    if logfilt:
        out[out<10**_np.ceil(_np.log10(out.mean()))] = 0
    return out
#eof
