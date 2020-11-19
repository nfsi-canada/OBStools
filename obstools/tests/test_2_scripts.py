import stdb
import numpy as np
from pkg_resources import resource_filename
from obspy.clients.fdsn import Client
from obstools.atacr import DayNoise, StaNoise, TFNoise
from obstools.atacr import EventStream, Power, Cross, Rotation
from obstools.atacr import utils, plotting

dbfile = resource_filename('obstools',
                           'examples/meta/M08A.pkl')

def test_download_data():
    from obstools.scripts import atacr_download_data as atacr
    args0 = atacr.get_daylong_arguments([
        dbfile, '--keys', '7D.M08A',
        '--start', '2012-03-08', '--end', '2012-03-10'])
    atacr.main(args=args0)
    atacr.main(args=args0)

def test_download_data_noP():
    from obstools.scripts import atacr_download_data as atacr
    args0 = atacr.get_daylong_arguments([
        dbfile, '--keys', '7D.M08A', '-O',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--channels', 'H'])
    atacr.main(args=args0)
    args0 = atacr.get_daylong_arguments([
        dbfile, '--keys', '7D.M08A',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--channels', 'H'])
    atacr.main(args=args0)

def test_download_data_noH():
    from obstools.scripts import atacr_download_data as atacr
    args0 = atacr.get_daylong_arguments([
        dbfile, '--keys', '7D.M08A', '-O',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--channels', 'P'])
    atacr.main(args=args0)
    args0 = atacr.get_daylong_arguments([
        dbfile, '--keys', '7D.M08A',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--channels', 'P'])
    atacr.main(args=args0)

def test_download_event():
    from obstools.scripts import atacr_download_event as atacr
    args0 = atacr.get_event_arguments([
        dbfile, '--keys', '7D.M08A',
        '--start', '2012-03-09', '--end', '2012-03-10',
        '--min-mag', '6.3', '--max-mag', '6.7'])
    atacr.main(args=args0)
    atacr.main(args=args0)

def test_download_event_noP():
    from obstools.scripts import atacr_download_event as atacr
    args0 = atacr.get_event_arguments([
        dbfile, '--keys', '7D.M08A', '-O',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--min-mag', '6.3', '--max-mag', '6.7',
        '--channels', 'H'])
    atacr.main(args=args0)
    args0 = atacr.get_event_arguments([
        dbfile, '--keys', '7D.M08A',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--min-mag', '6.3', '--max-mag', '6.7',
        '--channels', 'H'])
    atacr.main(args=args0)

def test_download_event_noH():
    from obstools.scripts import atacr_download_event as atacr
    args0 = atacr.get_event_arguments([
        dbfile, '--keys', '7D.M08A', '-O',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--min-mag', '6.3', '--max-mag', '6.7',
        '--channels', 'P'])
    atacr.main(args=args0)
    args0 = atacr.get_event_arguments([
        dbfile, '--keys', '7D.M08A',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--min-mag', '6.3', '--max-mag', '6.7',
        '--channels', 'P'])
    atacr.main(args=args0)

def test_dailyspec():
    from obstools.scripts import atacr_daily_spectra as atacr
    args0 = atacr.get_dailyspec_arguments([
        dbfile, '--keys', '7D.M08A', '-O', '--save-fig'])
    atacr.main(args=args0)

def test_cleanspec():
    from obstools.scripts import atacr_clean_spectra as atacr
    args0 = atacr.get_cleanspec_arguments([
        dbfile, '--keys', '7D.M08A', '-O', '--save-fig', '--figCross', '--figCoh'])
    atacr.main(args=args0)

def test_transfer():
    from obstools.scripts import atacr_transfer_functions as atacr
    args0 = atacr.get_transfer_arguments([
        dbfile, '--keys', '7D.M08A', '-O', '--save-fig', '--figTF'])
    atacr.main(args=args0)

def test_correct():
    from obstools.scripts import atacr_correct_event as atacr
    from obstools.scripts import atacr_download_event as event
    args0 = event.get_event_arguments([
        dbfile, '--keys', '7D.M08A', '-O',
        '--start', '2012-03-09', '--end', '2012-03-10',
        '--min-mag', '6.3', '--max-mag', '6.7'])
    event.main(args=args0)
    args0 = atacr.get_correct_arguments([
        dbfile, '--keys', '7D.M08A', '-O', '--figRaw', '--figClean', '--save-fig'])
    atacr.main(args=args0)

def test_comply_calculate():
    from obstools.scripts import comply_calculate as comply
    args0 = comply.get_comply_arguments([
        dbfile, '--keys', '7D.M08A', '-O', '--save-fig', '--fig'])
    comply.main(args=args0)

# def test_DayNoise():

#     trN1, trN2, trNZ, trNP = test_utils.test_get_data()
#     args = test_args.test_get_dailyspec_arguments()

#     for tr1, tr2, trZ, trP in zip(trN1, trN2, trNZ, trNP):

#         daynoise = DayNoise(tr1, tr2, trZ, trP,
#                             args.window, args.overlap, key='7D.M08A')

#         daynoise.QC_daily_spectra(pd=args.pd, tol=args.tol,
#                                   alpha=args.alpha, smooth=args.smooth)

#         daynoise.average_daily_spectra(
#             calc_rotation=True,
#             fig_average=True,
#             fig_coh_ph=True,
#             save='tmp', form='png')

#         plot = plotting(daynoise.f, daynoise.power,
#                         daynoise.power, daynoise.goodwins,
#                         daynoise.ncomp, key='7D.M08A')

#         daynoise.save('tmp')


# def test_StaNoise():

#     args = test_args.test_get_dailyspec_arguments()
#     stanoise = test_classes.test_stanoise_demo()

#     # Containers for power and cross spectra
#     coh_all = []
#     ph_all = []
#     coh_12_all = []
#     coh_1Z_all = []
#     coh_1P_all = []
#     coh_2Z_all = []
#     coh_2P_all = []
#     coh_ZP_all = []
#     ph_12_all = []
#     ph_1Z_all = []
#     ph_1P_all = []
#     ph_2Z_all = []
#     ph_2P_all = []
#     ph_ZP_all = []
#     ad_12_all = []
#     ad_1Z_all = []
#     ad_1P_all = []
#     ad_2Z_all = []
#     ad_2P_all = []
#     ad_ZP_all = []
#     nwins = []

#     for dn in stanoise.daylist:
#         dn.QC_daily_spectra()
#         dn.average_daily_spectra()
#         coh_all.append(dn.rotation.coh)
#         ph_all.append(dn.rotation.ph)

#         # Coherence
#         coh_12_all.append(
#             utils.smooth(
#                 utils.coherence(dn.cross.c12, dn.power.c11, dn.power.c22), 50))
#         coh_1Z_all.append(
#             utils.smooth(
#                 utils.coherence(dn.cross.c1Z, dn.power.c11, dn.power.cZZ), 50))
#         coh_1P_all.append(
#             utils.smooth(
#                 utils.coherence(dn.cross.c1P, dn.power.c11, dn.power.cPP), 50))
#         coh_2Z_all.append(
#             utils.smooth(
#                 utils.coherence(dn.cross.c2Z, dn.power.c22, dn.power.cZZ), 50))
#         coh_2P_all.append(
#             utils.smooth(
#                 utils.coherence(dn.cross.c2P, dn.power.c22, dn.power.cPP), 50))
#         coh_ZP_all.append(
#             utils.smooth(
#                 utils.coherence(dn.cross.cZP, dn.power.cZZ, dn.power.cPP), 50))

#         # Phase
#         try:
#             ph_12_all.append(
#                 180./np.pi*utils.phase(dn.cross.c12))
#         except:
#             ph_12_all.append(None)
#         try:
#             ph_1Z_all.append(
#                 180./np.pi*utils.phase(dn.cross.c1Z))
#         except:
#             ph_1Z_all.append(None)
#         try:
#             ph_1P_all.append(
#                 180./np.pi*utils.phase(dn.cross.c1P))
#         except:
#             ph_1P_all.append(None)
#         try:
#             ph_2Z_all.append(
#                 180./np.pi*utils.phase(dn.cross.c2Z))
#         except:
#             ph_2Z_all.append(None)
#         try:
#             ph_2P_all.append(
#                 180./np.pi*utils.phase(dn.cross.c2P))
#         except:
#             ph_2P_all.append(None)
#         try:
#             ph_ZP_all.append(
#                 180./np.pi*utils.phase(dn.cross.cZP))
#         except:
#             ph_ZP_all.append(None)

#         # Admittance
#         ad_12_all.append(utils.smooth(utils.admittance(
#             dn.cross.c12, dn.power.c11), 50))
#         ad_1Z_all.append(utils.smooth(utils.admittance(
#             dn.cross.c1Z, dn.power.c11), 50))
#         ad_1P_all.append(utils.smooth(utils.admittance(
#             dn.cross.c1P, dn.power.c11), 50))
#         ad_2Z_all.append(utils.smooth(utils.admittance(
#             dn.cross.c2Z, dn.power.c22), 50))
#         ad_2P_all.append(utils.smooth(utils.admittance(
#             dn.cross.c2P, dn.power.c22), 50))
#         ad_ZP_all.append(utils.smooth(utils.admittance(
#             dn.cross.cZP, dn.power.cZZ), 50))

#     # Convert to numpy arrays
#     coh_all = np.array(coh_all)
#     ph_all = np.array(ph_all)
#     coh_12_all = np.array(coh_12_all)
#     coh_1Z_all = np.array(coh_1Z_all)
#     coh_1P_all = np.array(coh_1P_all)
#     coh_2Z_all = np.array(coh_2Z_all)
#     coh_2P_all = np.array(coh_2P_all)
#     coh_ZP_all = np.array(coh_ZP_all)
#     ph_12_all = np.array(ph_12_all)
#     ph_1Z_all = np.array(ph_1Z_all)
#     ph_1P_all = np.array(ph_1P_all)
#     ph_2Z_all = np.array(ph_2Z_all)
#     ph_2P_all = np.array(ph_2P_all)
#     ph_ZP_all = np.array(ph_ZP_all)
#     ad_12_all = np.array(ad_12_all)
#     ad_1Z_all = np.array(ad_1Z_all)
#     ad_1P_all = np.array(ad_1P_all)
#     ad_2Z_all = np.array(ad_2Z_all)
#     ad_2P_all = np.array(ad_2P_all)
#     ad_ZP_all = np.array(ad_ZP_all)

#     # Store transfer functions as objects for plotting
#     coh = Cross(coh_12_all, coh_1Z_all, coh_1P_all,
#                 coh_2Z_all, coh_2P_all, coh_ZP_all)
#     ph = Cross(ph_12_all, ph_1Z_all, ph_1P_all,
#                ph_2Z_all, ph_2P_all, ph_ZP_all)
#     ad = Cross(ad_12_all, ad_1Z_all, ad_1P_all,
#                ad_2Z_all, ad_2P_all, ad_ZP_all)

#     stanoise.QC_sta_spectra(pd=args.pd, tol=args.tol, alpha=args.alpha,
#                             fig_QC=True, debug=False)

#     stanoise.average_sta_spectra(fig_average=True)

#     plot = plotting.fig_av_cross(stanoise.f, coh,
#                                  stanoise.gooddays,
#                                  'Coherence', stanoise.ncomp,
#                                  key='7D.M08A', lw=0.5)
#     plot.close()

#     plot = plotting.fig_coh_ph(coh_all, ph_all, stanoise.direc)
#     plot.close()

#     return stanoise


# def test_TFNoise(tmp_path):

#     args = test_args.test_get_correct_arguments()
#     tfnoise_day = test_classes.test_tfnoise_day_demo(tmp_path)
#     tfnoise_sta = test_classes.test_tfnoise_sta_demo(tmp_path)
#     f = tfnoise_day.f

#     plot = plotting.fig_TF(
#         f, [tfnoise_day.transfunc],
#         tfnoise_day.tf_list,
#         tfnoise_sta.transfunc,
#         tfnoise_sta.tf_list,
#         skey='7D.M08A')
#     plot.close()

#     eventstream = test_classes.test_evstream_demo(tmp_path)

#     plot = plotting.fig_event_raw(
#         eventstream, fmin=args.fmin, fmax=args.fmax)
#     plot.close()

#     eventstream.correct_data(tfnoise_day)
#     correct = eventstream.correct
#     plot = plotting.fig_event_corrected(
#         eventstream, tfnoise_day.tf_list)
#     plot.close()

#     eventstream = test_classes.test_evstream_demo(tmp_path)
#     eventstream.correct_data(tfnoise_sta)
#     plot = plotting.fig_event_corrected(
#         eventstream, tfnoise_sta.tf_list)
#     plot.close()

