'''
MODULE obs_proc.py

Main processing functions for tilt and compliance correction
of vertical component OBS data. Algorithm is based on

Bell, S.W., Forsyth, D.W., and Ruan, Y., Removing noise from
the vertical component records of ocean-bottom seismometers: Results
from Year one of the Cascadia Initiative, Bull. Seism. Soc. Am., 
105, 300-313, 2015.

---
Pascal Audet
pascal.audet@uottawa

Last updated: 17 November 2015

'''

import numpy as np
from matplotlib import pyplot as plt
from obs import obs_utils as obsu
from obs import obs_plot as obspl

def obs_proc_tilt(trN1, trN2, trNZ, tr1, tr2, trZ, plot=True):
    """
    Function to process OBS data for tilt correction
    by applying a transfer function from the noise records.

    Tilt calculation must be performed prior to this step. Tilt
    value is stored in the header of the noise trace trNZ.
    See obs_util.py for details.

    :return trNZ_p, trZ_p: tilt-free vertical traces as traces

    """

    # Get parameters from traces
    nN = trNZ.stats.npts
    dN = trNZ.stats.delta
    nn = trZ.stats.npts
    dn = trZ.stats.delta
    tilt = trNZ.stats.sac.user7
    ws = int(trNZ.stats.sac.user8)
    ss = int(trNZ.stats.sac.user9)

    # Frequency bounds for transfer function
    f1 = 0.005
    f2 = 0.06

    # Time and frequency axes
    timeN = np.arange(nN)*dN
    time = np.arange(nn)*dn
    freq = np.fft.fftfreq(ws, dN)

    # Copy traces
    tr_N1 = trN1.copy()
    tr_1 = tr1.copy()

    # Rotate H1 to tilt direction
    tr_N1.data = obsu.rotate_dir(trN1.data, trN2.data, tilt)
    tr_1.data = obsu.rotate_dir(tr1.data, tr2.data, tilt)

    if plot:
        obspl.obs_plot_2traces(trNZ, tr_N1, title='Noise traces Z, rotated H1')
    
    # Get transfer functions between 1 and Z
    Ad, Co, Ph, eAd, eCo, ePh = obsu.calculate_trfs(tr_N1, trNZ)

    if plot:
        obspl.obs_plot_trf_tilt(Ad, Co, Ph, eAd, eCo, ePh, freq)

    # Get analytical transfer functions
    pAd, pPh = obsu.get_trfs_lsqfit(Ad, Co, Ph, eAd, eCo, ePh, freq, f1, f2, plot=plot)

    # Get predicted Z from transfer functions
    trZ_p = obsu.get_trZ_p(tr_1, pAd, pPh, f1, f2)
    trNZ_p = obsu.get_trZ_p(tr_N1, pAd, pPh, f1, f2)

    # Remove predicted Z from observed
    trZ_pp = trZ.copy()
    trZ_pp.data = trZ.data - trZ_p.data

    trNZ_pp = trNZ.copy()
    trNZ_pp.data = trNZ.data - trNZ_p.data

    if plot:
        obspl.obs_plot_3traces(trZ, trZ_p, trZ_pp, \
                title='Observed, predicted and corrected vertical seismogram')
        
    return trNZ_pp, trZ_pp


def obs_proc_compliance(trNZ, trNP, trZ, trP, plot=True):
    """
    Function to process OBS data for compliance correction
    by applying a transfer function from the noise records.

    Tilt calculation must be performed prior to this step. Tilt
    value is stored in the header of the noise trace trNZ.
    See obs_util.py for details.

    :return trNZ_p, trZ_p: tilt-free vertical traces as traces

    """

    # Get parameters from traces
    nN = trNZ.stats.npts
    dN = trNZ.stats.delta
    nn = trZ.stats.npts
    dn = trZ.stats.delta
    ws = int(trNZ.stats.sac.user8)
    ss = int(trNZ.stats.sac.user9)
    trNP.data *= -1.
    trP.data *= -1.

    # Frequency bounds for transfer function
    f1 = 0.005

    # Second frequency bound based on station depth
    d = trNZ.stats.sac.stel
    f2 = np.sqrt(9.81/2./np.pi/d)
    print('Max frequency for compliance correction')
    print('based on station depth:',f2,' Hz')

    # Time and frequency axes
    timeN = np.arange(nN)*dN
    time = np.arange(nn)*dn
    freq = np.fft.fftfreq(ws, dN)

    if plot:
        obspl.obs_plot_2traces(trNZ, trNP, title='Noise traces Z, P')

    
    # Get transfer functions between P and Z
    Ad, Co, Ph, eAd, eCo, ePh = obsu.calculate_trfs(trNP, trNZ)

    if plot:
        obspl.obs_plot_trf_compliance(Ad, Co, Ph, eAd, eCo, ePh, freq)

    # Get analytical transfer functions
    pAd, pPh = obsu.get_trfs_lsqfit(Ad, Co, Ph, eAd, eCo, ePh, freq, f1, f2, plot=plot)

    # Get predicted Z from transfer functions
    trZ_p = obsu.get_trZ_p(trP, pAd, pPh, f1, f2)

    # Remove predicted Z from observed
    trZ_pp = trZ.copy()
    trZ_pp.data = trZ.data - trZ_p.data/2.

    if plot:
        obspl.obs_plot_3traces(trZ, trZ_p, trZ_pp, f1=0.1, f2=1.0, \
                title='Observed, predicted and corrected vertical seismogram')

    return trZ_pp 


