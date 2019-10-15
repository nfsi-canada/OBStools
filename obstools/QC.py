from scipy.signal import spectrogram, savgol_filter, detrend
from scipy.linalg import norm
from obstools import utils
import matplotlib.pyplot as plt
import numpy as np


def QC_daily_spectra(tr1, tr2, trZ, trP, window, overlap, 
    pd=[0.004, 2.0], tol=1.5, alpha=0.05, smooth=True, fig_QC_specs=False, debug=False):

    # Points in window
    ws = int(window/tr1.stats.delta)

    # Number of points to overlap
    ss = int(window*overlap)/tr1.stats.delta

    # Sampling frequency
    fs = tr1.stats.sampling_rate

    # Get spectrograms for single day-long keys
    f, t, psd1 = spectrogram(tr1.data, fs, window='hanning', nperseg=ws, noverlap=ss)
    f, t, psd2 = spectrogram(tr2.data, fs, window='hanning', nperseg=ws, noverlap=ss)
    f, t, psdZ = spectrogram(trZ.data, fs, window='hanning', nperseg=ws, noverlap=ss)
    f, t, psdP = spectrogram(trP.data, fs, window='hanning', nperseg=ws, noverlap=ss)

    if debug:
        plt.figure(1)
        plt.subplot(4,1,1)
        plt.pcolormesh(t, f, np.log(psd1))
        plt.subplot(4,1,2)
        plt.pcolormesh(t, f, np.log(psd2))
        plt.subplot(4,1,3)
        plt.pcolormesh(t, f, np.log(psdZ))
        plt.subplot(4,1,4)
        plt.pcolormesh(t, f, np.log(psdP))
        plt.tight_layout()
        plt.show()

    # Select bandpass frequencies
    f_wt = f[(f>pd[0]) & (f<pd[1])]

    if smooth:
        # Smooth out the log of the PSDs
        sl_psd1 = savgol_filter(np.log(psd1[(f>pd[0]) & (f<pd[1]),:]), 51, 3, axis=0)
        sl_psd2 = savgol_filter(np.log(psd2[(f>pd[0]) & (f<pd[1]),:]), 51, 3, axis=0)
        sl_psdZ = savgol_filter(np.log(psdZ[(f>pd[0]) & (f<pd[1]),:]), 51, 3, axis=0)
        sl_psdP = savgol_filter(np.log(psdP[(f>pd[0]) & (f<pd[1]),:]), 51, 3, axis=0)
    else:
        # Take the log of the PSDs
        sl_psd1 = np.log(psd1[(f>pd[0]) & (f<pd[1]),:])
        sl_psd2 = np.log(psd2[(f>pd[0]) & (f<pd[1]),:])
        sl_psdZ = np.log(psdZ[(f>pd[0]) & (f<pd[1]),:])
        sl_psdP = np.log(psdP[(f>pd[0]) & (f<pd[1]),:])

    # Remove mean of the log PSDs
    dsl_psd1 = sl_psd1 - np.mean(sl_psd1, axis=0)
    dsl_psd2 = sl_psd2 - np.mean(sl_psd2, axis=0)
    dsl_psdZ = sl_psdZ - np.mean(sl_psdZ, axis=0)
    dsl_psdP = sl_psdP - np.mean(sl_psdP, axis=0)

    if debug:
        plt.figure(2)
        plt.subplot(4,1,1)
        plt.semilogx(f_wt, sl_psd1, 'r', lw=0.5)
        plt.subplot(4,1,2)
        plt.semilogx(f_wt, sl_psd2, 'b', lw=0.5)
        plt.subplot(4,1,3)
        plt.semilogx(f_wt, sl_psdZ, 'g', lw=0.5)
        plt.subplot(4,1,4)
        plt.semilogx(f_wt, sl_psdP, 'k', lw=0.5)
        plt.tight_layout()
        plt.show()

    # Cycle through to kill high-std-norm windows
    moveon = False
    goodwins = np.repeat([True], len(t))

    while moveon == False:
        ubernorm = np.empty((4, np.sum(goodwins)))
        for ind_u, dsl in enumerate([dsl_psd1, dsl_psd2, dsl_psdZ, dsl_psdP]):
            normvar = np.zeros(np.sum(goodwins))
            for ii in range(0, np.sum(goodwins)):
                ind = np.argwhere(goodwins); ind[ii] = False
                normvar[ii] = norm(np.std(dsl[:, ind], axis=1), ord=2)
            ubernorm[ind_u, :] = np.median(normvar) - normvar 

        penalty = np.sum(detrend(ubernorm, type='constant'), axis=0)

        if debug:
            plt.figure(4)
            plt.plot(range(0,np.sum(goodwins)), detrend(ubernorm, type='constant')[0], 'o-')
            plt.plot(range(0,np.sum(goodwins)), detrend(ubernorm, type='constant')[1], 'o-')
            plt.plot(range(0,np.sum(goodwins)), detrend(ubernorm, type='constant')[2], 'o-')
            plt.plot(range(0,np.sum(goodwins)), detrend(ubernorm, type='constant')[3], 'o-')
            plt.show()
            plt.figure(5)
            plt.plot(range(0,np.sum(goodwins)), np.sum(ubernorm, axis=0), 'o-')
            plt.show()

        kill = penalty > tol*np.std(penalty)
        if np.sum(kill)==0: return
        trypenalty = penalty[np.argwhere(kill == False)].T[0]

        if utils.ftest(penalty, 1, trypenalty, 1) < alpha:
            goodwins[np.argwhere(kill==True)] = False
            moveon = False
        else:
            moveon = True

    if fig_QC_specs:
        plt.figure(6)
        plt.subplot(4,1,1)
        plt.semilogx(f_wt, sl_psd1[:,goodwins], 'k', lw=0.5)
        if np.sum(~goodwins)>0:
            plt.semilogx(f_wt, sl_psd1[:,~goodwins], 'r', lw=0.5)
        plt.title('H1 windows', fontdict={'fontsize': 8})
        plt.subplot(4,1,2)
        plt.semilogx(f_wt, sl_psd2[:,goodwins], 'k', lw=0.5)
        if np.sum(~goodwins)>0:
            plt.semilogx(f_wt, sl_psd2[:,~goodwins], 'r', lw=0.5)
        plt.title('H2 windows', fontdict={'fontsize': 8})
        plt.subplot(4,1,3)
        plt.semilogx(f_wt, sl_psdZ[:,goodwins], 'k', lw=0.5)
        if np.sum(~goodwins)>0:
            plt.semilogx(f_wt, sl_psdZ[:,~goodwins], 'r', lw=0.5, label='')
        plt.title('HZ windows', fontdict={'fontsize': 8})
        plt.subplot(4,1,4)
        plt.semilogx(f_wt, sl_psdP[:,goodwins], 'k', lw=0.5)
        if np.sum(~goodwins)>0:
            plt.semilogx(f_wt, sl_psdP[:,~goodwins], 'r', lw=0.5)
        plt.title('HP windows', fontdict={'fontsize': 8})
        plt.tight_layout()
        plt.show()

    return goodwins


def spec_calc_daily_spectra(tr1, tr2, trZ, trP, goodwins, window, overlap, \
    calc_rotation=True, fig_av_specs=False, debug=False):

    # Points in window
    ws = int(window/tr1.stats.delta)

    # Number of points in step
    ss = int(window*(1.-overlap)/tr1.stats.delta)

    ft1, f = utils.calculate_windowed_fft(tr1, ws, ss)
    ft2, f = utils.calculate_windowed_fft(tr2, ws, ss)
    ftZ, f = utils.calculate_windowed_fft(trZ, ws, ss)
    ftP, f = utils.calculate_windowed_fft(trP, ws, ss)

    cc1 = np.abs(np.mean(ft1[goodwins,:]*np.conj(ft1[goodwins,:]), axis=0))
    cc2 = np.abs(np.mean(ft2[goodwins,:]*np.conj(ft2[goodwins,:]), axis=0))
    ccZ = np.abs(np.mean(ftZ[goodwins,:]*np.conj(ftZ[goodwins,:]), axis=0))
    ccP = np.abs(np.mean(ftP[goodwins,:]*np.conj(ftP[goodwins,:]), axis=0))
    if np.sum(~goodwins) > 0:
        bcc1 = np.abs(np.mean(ft1[~goodwins,:]*np.conj(ft1[~goodwins,:]), axis=0))
        bcc2 = np.abs(np.mean(ft2[~goodwins,:]*np.conj(ft2[~goodwins,:]), axis=0))
        bccZ = np.abs(np.mean(ftZ[~goodwins,:]*np.conj(ftZ[~goodwins,:]), axis=0))
        bccP = np.abs(np.mean(ftP[~goodwins,:]*np.conj(ftP[~goodwins,:]), axis=0))

    c12 = np.mean(ft1[goodwins,:]*np.conj(ft2[goodwins,:]), axis=0)
    c1Z = np.mean(ft1[goodwins,:]*np.conj(ftZ[goodwins,:]), axis=0)
    c1P = np.mean(ft1[goodwins,:]*np.conj(ftP[goodwins,:]), axis=0)
    c2Z = np.mean(ft2[goodwins,:]*np.conj(ftZ[goodwins,:]), axis=0)
    c2P = np.mean(ft2[goodwins,:]*np.conj(ftP[goodwins,:]), axis=0)
    cZP = np.mean(ftZ[goodwins,:]*np.conj(ftP[goodwins,:]), axis=0)

    if fig_av_specs:
        plt.figure()
        plt.subplot(411)
        plt.semilogx(f, savgol_filter(np.log(cc1[0:len(f)]), 51, 3), 'k', lw=0.5)
        if np.sum(~goodwins)>0:
            plt.semilogx(f, savgol_filter(np.log(bcc1[0:len(f)]), 51, 3), 'r', lw=0.5)
        plt.title('Daily average H1', fontdict={'fontsize': 8})
        plt.subplot(412)
        plt.semilogx(f, savgol_filter(np.log(cc2[0:len(f)]), 51, 3), 'k', lw=0.5)
        if np.sum(~goodwins)>0:
            plt.semilogx(f, savgol_filter(np.log(bcc2[0:len(f)]), 51, 3), 'r', lw=0.5)
        plt.title('Daily average H2', fontdict={'fontsize': 8})
        plt.subplot(413)
        plt.semilogx(f, savgol_filter(np.log(ccZ[0:len(f)]), 51, 3), 'k', lw=0.5)
        if np.sum(~goodwins)>0:
            plt.semilogx(f, savgol_filter(np.log(bccZ[0:len(f)]), 51, 3), 'r', lw=0.5)
        plt.title('Daily average HZ', fontdict={'fontsize': 8})
        plt.subplot(414)
        plt.semilogx(f, savgol_filter(np.log(ccP[0:len(f)]), 51, 3), 'k', lw=0.5)
        if np.sum(~goodwins)>0:
            plt.semilogx(f, savgol_filter(np.log(bccP[0:len(f)]), 51, 3), 'r', lw=0.5)
        plt.title('Daily average HP', fontdict={'fontsize': 8})
        plt.tight_layout()
        plt.show()

    if calc_rotation:
        coh, ph, tilt, coh_value, phase_value = utils.calculate_tilt( \
            ft1, ft2, ftZ, f, goodwins, plot=False)


