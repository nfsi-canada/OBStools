import numpy as _np
import obspy
from scipy.signal import stft, detrend
import copy
import obstools
import matplotlib.pyplot as plt

class OBSMetrics(object):

        def __init__(self,tr1=None, tr2=None, trZ=None, trP=None,overlap=0.3,csd=None,f=None):
                self.overlap = overlap
                self.traces=dict()
                self.csd_pairs = ['1P','2P','1Z','2Z','ZP','PP','11','22','ZZ']
                self.traces['1'],self.traces['2'],self.traces['Z'],self.traces['P'] = [],[],[],[]
                self.csd = dict()
                self.csd['A'] = dict()
                self.csd['B'] = dict()
                self.csd['AB'] = dict()
                for p in self.csd_pairs:
                        self.csd['A'][p] = []
                        self.csd['B'][p] = []
                        self.csd['AB'][p] = []
                self.ntraces = 0
                self.StaNoise = None
                self.f = f
                if tr1 is not None:
                        self._add_traces(tr1=tr1,tr2=tr2,trZ=trZ,trP=trP)
                        self.A = self.traces.copy()
                        self.B = self.traces.copy()
                        self._meta()
                        self.updatespec()
                if csd is not None:
                        if isinstance(csd,obstools.atacr.StaNoise):
                                dct = csd.power.__dict__.copy()
                                dct.update(csd.cross.__dict__.copy())
                                self.StaNoise = csd
                                csd = dct.copy()
                        for a in self.csd_pairs:
                                try:
                                        self.csd['A'][a] = csd['c'+a]
                                        self.csd['B'][a] = csd['c'+a]
                                        self.csd['AB'][a] = csd['c'+a]
                                except:
                                        try:
                                                self.csd['A'][a] = csd['c'+a[::-1]]
                                                self.csd['B'][a] = csd['c'+a[::-1]]
                                                self.csd['AB'][a] = csd['c'+a[::-1]]
                                        except:
                                                try:
                                                        self.csd['A'][a] = csd[a]
                                                        self.csd['B'][a] = csd[a]
                                                        self.csd['AB'][a] = csd[a]
                                                except:
                                                        self.csd['A'][a] = csd[a[::-1]]
                                                        self.csd['B'][a] = csd[a[::-1]]
                                                        self.csd['AB'][a] = csd[a[::-1]]
        def ft(self,r,window=None,overlap=None):
                tr = self.traces[r]
                ft = []
                for tr in self.traces[r]:
                        f,_t,ft0 = self._stft(tr,window=window,overlap=overlap)
                        ft.append(_np.mean(ft0,axis=-1))
                return f,_np.array(ft).squeeze()
        def psd(self,r,window=None,overlap=None):
                tr = self.traces[r]
                ft = []
                f = self.f
                for tr in self.traces[r]:
                        f,ft0 = self._psd(tr,window=window,overlap=overlap)
                        ft.append(ft0)
                ft = _np.mean(_np.array(ft).squeeze(),axis=0)
                ft = ft[f>=0]
                f = f[f>=0]
                return f,ft
        def Phase(self,r,s=False):
                '''
                Takes a single NDarray containing the cross-psd (complex) between two components and pulls the phase out of the imaginary component of the array. Output is in degrees for a quadrature circle (-180,180).
                '''
                # r = ''.join(sorted(r))
                f = self.f[self.f>=0]
                ph = [self.calc_phase(ab) for ab in self.csd['AB'][r]]
                ph = _np.real(_np.array(ph).squeeze()[self.f>=0])
                if s:
                        ph = self.smooth(ph)
                ph = ph[f>=0]
                f = f[f>=0]
                return f,ph
        def Admittance(self,r,s=False):
                '''
                Takes two NDarrays of equal shape, the first (a) containing the cross-power-spectral density between the two components and the second (b) containing the auto-power-spectral density of the secondary component.
                ie, For the spectral admittance between Z and P, a is the cross psd between the two and b is the auto-psd of P.
                '''
                # r = ''.join(sorted(r))
                f = self.f[self.f>=0]
                adm = [self.calc_admittance(self.csd['AB'][r][i],self.csd['B'][r[1] + r[1]][i]) for i in range(len(self.csd['AB'][r]))]
                adm = _np.real(_np.array(adm).squeeze()[self.f>=0])
                if s:
                        adm = self.smooth(adm)
                adm = adm[f>=0]
                f = f[f>=0]
                return f,adm
        def Coherence(self,r,s=False):
                '''
                Takes three NDarrays of equal shape, the first (ab) containing the cross-power-spectral density between the two components, the second and third (aa and bb) containing the auto-power-spectral density of these two components.
                ie, For the spectral coherence between Z and P, a is the cross psd between the two and b and c are the auto-psd of Z and P, respectively.
                '''
                # r = ''.join(sorted(r))
                f = self.f[self.f>=0]
                coh = [self.calc_coherence(   self.csd['AB'][r][i],    self.csd['A'][r[0] + r[0]][i],   self.csd['B'][r[1] + r[1]][i])     for i in range(len(self.csd['AB'][r]))]
                coh = _np.abs(_np.array(coh).squeeze()[self.f>=0])
                if s:
                        # coh = _np.array([self.smooth(c) for c in coh])
                        coh = self.smooth(coh)
                coh = coh[f>=0]
                f = f[f>=0]
                return f,coh
        def CrossSpec(self,A,B,fs=None,return_onesided=False,window=None):
                if isinstance(A,obspy.core.trace.Trace):
                        fs = 1/A.stats.delta
                        A = A.data
                if isinstance(B,obspy.core.trace.Trace):
                        fs = 1/B.stats.delta
                        B = B.data
                f,spec_AB= self._csd_helper(A,B,window=window)
                f,spec_AA= self._csd_helper(A,A,window=window)
                f,spec_BB= self._csd_helper(B,B,window=window)
                self.f = f
                return f,spec_AB,_np.abs(spec_AA),_np.abs(spec_BB)
        def Metrics(self,r):
                _,coh = self.Coherence(r)
                _,adm = self.Admittance(r)
                _,ph = self.Phase(r)
                return self.f,coh,adm,ph
        def smooth(self,d,k=10):
                return _np.convolve(d, _np.ones(k) / k, mode='same')
        def updatespec(self,window=None,overlap=None):
                self.csd = dict()
                self.csd['A'] = dict()
                self.csd['B'] = dict()
                self.csd['AB'] = dict()
                for p in self.csd_pairs:
                        self.csd['A'][p] = []
                        self.csd['B'][p] = []
                        self.csd['AB'][p] = []
                for i in range(len(self.traces['Z'])):
                        for p in self.csd_pairs:
                                A = self.A[p[0]][i].copy()
                                B = self.B[p[1]][i].copy()
                                f,spec_AB= self._csd_helper(A.data,B.data,window=window,overlap=overlap)
                                self.csd['AB'][p].append(spec_AB)
                                A = self.A[p[0]][i].copy()
                                B = self.A[p[1]][i].copy()
                                f,spec_AB= self._csd_helper(A.data,B.data,window=window,overlap=overlap)
                                self.csd['A'][p].append(spec_AB)
                                A = self.B[p[0]][i].copy()
                                B = self.B[p[1]][i].copy()
                                f,spec_AB= self._csd_helper(A.data,B.data,window=window,overlap=overlap)
                                self.csd['B'][p].append(spec_AB)
                self.f = f
        def spectrogram(self,r,i=0,window=None,overlap=None,plot=False,yscale='linear',cmap='magma'):
                f,t,s = self._stft(self.traces[r][i].data,scaling='psd',window=window,overlap=overlap,return_onesided=True)
                f,t,s = f.T,t.T,s.T
                s = _np.abs(s)
                if plot:
                        fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(17,10),sharex='all',layout='constrained')
                        x = self.traces['Z'][0].times()
                        y = self.traces['Z'][0].data
                        ax[0].plot(x,y)
                        ax[0].set_xlim(x[0],x[-1])
                        ax[1].set_ylabel('Frequency (Hz)',fontweight='bold')
                        ax[1].set_xlabel('Time (s)',fontweight='bold')
                        pc = ax[1].pcolormesh(t,f, 10*_np.log10(s), cmap = cmap, shading= 'auto')
                        ax[1].set_ylabel('Frequency (Hz)',fontweight='bold')
                        ax[1].set_xlabel('Time (s)',fontweight='bold')
                        ax[1].set_yscale(yscale)
                        ax[1].set_ylim(f[1],f[-1])
                        [ax[rw].axvline(t[0],color='k') for rw in range(2)],[ax[rw].axvline(t[-1],color='k') for rw in range(2)]
                        ax[0].set_xlim(t[0],t[-1])
                        fig.colorbar(pc, ax=ax[1], pad=0.01, label='dB')
                        ttl = self.traces[r][i].stats.station +  ' | ' + self.traces[r][i].stats.location
                        fig.suptitle(ttl,fontweight='bold')
                return f,t,s
        def _csd_helper(self,a,b,window=None,overlap=None):
                f,_t,a_ft = self._stft(a,window=window,overlap=overlap)
                f,_t,b_ft = self._stft(b,window=window,overlap=overlap)
                cab = _np.mean(self._csd(a_ft,b_ft),axis=0)
                return f,cab
        def _stft(self,tr,scaling='spectrum',window=None,overlap=None,return_onesided=False):
                wind,ws,ss = self._window(window=window,overlap=overlap)
                _f, _t, ft = stft(tr,self.fs, return_onesided=return_onesided, boundary=None,padded=False, window=wind, nperseg=ws, noverlap=ss,detrend='constant',scaling=scaling)
                _f = _f.reshape(-1)
                ft = _np.atleast_2d(ft)
                ft = ft*ws
                # ft = _np.mean(ft,axis=-1)
                return _f.T,_t.T,ft.T
        def _window(self,window=None,overlap=None):
                if overlap is None:
                        overlap=self.overlap
                if window is None:
                        window=self.window
                # Points in window
                ws = int(window/self.dt)
                # Number of points to overlap
                ss = int(window*overlap/self.dt)
                # hanning window
                hanning = _np.hanning(2*ss)
                wind = _np.ones(ws)
                wind[0:ss] = hanning[0:ss]
                wind[-ss:ws] = hanning[ss:ws]
                return wind,ws,ss
        def _csd(self,a_ft,b_ft):
                cab = a_ft*_np.conj(b_ft)
                return cab
        def _add_traces(self,tr1=None,tr2=None,trZ=None,trP=None):
                self.traces['1'].append(tr1.copy())
                self.traces['2'].append(tr2.copy())
                self.traces['Z'].append(trZ.copy())
                self.traces['P'].append(trP.copy())
                self.A,self.B = self.traces.copy(),self.traces.copy()
                self.ntraces+=1
        def _psd(self,tr,window=None,overlap=None):
                f,_t,ft = self._stft(tr,scaling='psd',return_onesided=True,window=window,overlap=overlap)
                psd = _np.abs(ft)**2
                return f,psd
        def calc_phase(self,ab):
                ph = _np.angle(ab,deg=True)
                return ph
        def calc_admittance(self,ab,bb):
                ad = _np.abs(ab)/bb
                return ad
        def calc_coherence(self,ab,aa,bb):
                coh = ((abs(ab)**2)/(aa*bb))
                return coh
        def __add__(self,other):
                if isinstance(other,OBSMetrics):
                        self._add_traces(tr1=other.traces['1'][0],tr2=other.traces['2'][0],trZ=other.traces['Z'][0],trP=other.traces['P'][0])
                elif isinstance(other,Stream):
                        self._add_traces(tr1=other.select(component='1'),tr2=other.select(component='2'),trZ=other.select(component='Z'),trP=other.select(component='H'))
                self._meta()
                self.updatespec()
                return self
        def copy(self):
                return copy.deepcopy(self)
        def append(self,other):
                self.__add__(other)
        def __truediv__(self, other):
                self = self.copy()
                if isinstance(other,OBSMetrics):
                        self.B['1']=other.traces['1'].copy()
                        self.B['2']=other.traces['2'].copy()
                        self.B['Z']=other.traces['Z'].copy()
                        self.B['P']=other.traces['P'].copy()
                elif isinstance(other,Stream):
                        self.B['1'] = other.select(component='1').copy()
                        self.B['2'] = other.select(component='2').copy()
                        self.B['Z'] = other.select(component='Z').copy()
                        self.B['P'] = other.select(component='H').copy()
                self.updatespec()
                return self.copy()
        def _meta(self):
                self.dt = _np.unique([A.stats.delta for A in [self.traces['1'][0],self.traces['2'][0],self.traces['Z'][0],self.traces['P'][0]] if isinstance(A,obspy.core.trace.Trace)])
                self.fs = 1/self.dt
                self.tlen = _np.size(self.traces['Z'][0])/self.fs
                if self.tlen>7200:
                        self.window = 7200
                else:
                        self.window = 500