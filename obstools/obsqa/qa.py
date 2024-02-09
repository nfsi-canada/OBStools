from ObsQA.imports import *
def _next_pow2(n):
        return int(round(2**np.ceil(np.log2(n))))

def Coherence(a,b,c):
        '''
        Takes three NDarrays of equal shape, the first (a) containing the cross-power-spectral density between the two components, the second and third (b and c) containing the auto-power-spectral density of the two components.
        ie, For the spectral coherence between Z and P, a is the cross psd between the two and b and c are the auto-psd of Z and P, respectively.
        '''
        if not (a.shape==b.shape==c.shape):
                raise Exception("Input arrays are not the same shapes")
        coh = ((abs(a)**2)/(b*c))
        return coh

def Phase(a):
        '''
        Takes a single NDarray containing the cross-psd (complex) between two components and pulls the phase out of the imaginary component of the array. Output is in degrees for a quadrature circle (-180,180).
        '''
        ph = 180/np.pi*np.arctan2(np.imag(a),np.real(a))
        return ph

def Admittance(a,b):
        '''
        Takes two NDarrays of equal shape, the first (a) containing the cross-power-spectral density between the two components and the second (b) containing the auto-power-spectral density of the secondary component.
        ie, For the spectral admittance between Z and P, a is the cross psd between the two and b is the auto-psd of P.
        '''
        if not (a.shape==b.shape):
                raise Exception("Input arrays are not the same shapes")
        ad = np.abs(a)/b
        return ad

def PZCOH(P,Z):
        fs = 1/P.stats.delta
        nfft = _next_pow2(len(Z.data))
        f,spec_zp= csd(Z.data,P.data, fs=fs,return_onesided=True,nfft=nfft)
        f,spec_zz= csd(Z.data,Z.data, fs=fs,return_onesided=True,nfft=nfft)
        f,spec_pp= csd(P.data,P.data, fs=fs,return_onesided=True,nfft=nfft)
        coh = Coherence(spec_zp,spec_zz,spec_pp)
        return f,coh

#eof