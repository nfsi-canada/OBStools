from ObsQA.imports import *

def DownloadEvents(catalog,netsta_names=None,Minmag=6.3,Maxmag=6.7,limit=1000,pre_event_min_aperture=1,subfolder=None):
        # sys.stdout.flush()
        if subfolder is not None:
                logoutput = subfolder + '_Step_2_7_EventDownload_logfile.log'
        else:
                logoutput = '_Step_2_7_EventDownload_logfile.log'
        datafolder = './Data/'
        staquery_output = datafolder + 'sta_query.pkl'
        dateformat = '%Y.%j.%H.%M'
        print('----Begin Event Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i in range(len(catalog)):
                csta = catalog.iloc[i]
                S = csta['Station']#station
                N = csta['Network']
                staname = str(N) + '.' + str(S)
                if os.path.isdir(datafolder + staname)==False:
                        os.system('mkdir ' + datafolder + '/' + staname)
                ObsQA.io.build_staquery(catalog[(catalog.Network==N) & (catalog.Station==S)],staquery_output)
                log_fout = datafolder + staname + '/' + logoutput
                original = sys.stdout
                sys.stdout = open(log_fout,'w+')
                print('--' + staname + '--',flush=True)
                for j in range(len(csta.Events)):
                        print(staname + ' Station ' +str(i+1) + '/' + str(len(catalog)) + ' - Event ' + str(j+1) + '/' + str(len(csta.Events)),flush=True)
                        ev = csta.Events[j]
                        EventStart = UTCDateTime.strptime(ev,dateformat)
                        EventEnd = UTCDateTime.strptime(ev,dateformat) + datetime.timedelta(minutes=pre_event_min_aperture)
                        args = [staquery_output,'--start={}'.format(EventStart), '--end={}'.format(EventEnd),'--min-mag={}'.format(Minmag),'--max-mag={}'.format(Maxmag),'--limit={}'.format(limit)]
                        # with open(log_fout, 'w') as sys.stdout:
                        atacr_download_event.main(atacr_download_data.get_event_arguments(args))
        print(' ')
        print('----Event Download Complete----')
        sys.stdout = original

def DownloadNoise(catalog,netsta_names=None,pre_event_day_aperture=30,subfolder=None):
        # sys.stdout.flush()
        if subfolder is not None:
                logoutput = subfolder + '_Step_3_7_NoiseDownload_logfile.log'
        else:
                logoutput = '_Step_3_7_NoiseDownload_logfile.log'
        datafolder = './Data/'
        staquery_output = datafolder + 'sta_query.pkl'
        dateformat = '%Y.%j.%H.%M'

        print('----Begin Noise Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i in range(len(catalog)):
                csta = catalog.iloc[i]
                S = csta['Station']#station
                N = csta['Network']
                staname = str(N) + '.' + str(S)
                if os.path.isdir(datafolder + staname)==False:
                        os.system('mkdir ' + datafolder + '/' + staname)
                ObsQA.io.build_staquery(catalog[(catalog.Network==N) & (catalog.Station==S)],staquery_output)
                log_fout = datafolder + staname + '/' + logoutput
                # with open(log_fout, 'w') as sys.stdout:
                original = sys.stdout
                sys.stdout = open(log_fout,'w+')
                print('--' + staname + '--',flush=True)
                for j in range(len(csta.Events)):
                        print(staname + ' Station ' +str(i+1) + '/' + str(len(catalog)) + ' - Event ' + str(j+1) + '/' + str(len(csta.Events)),flush=True)
                        ev = csta.Events[j]
                        NoiseStart = UTCDateTime.strptime(ev,dateformat) - datetime.timedelta(days=pre_event_day_aperture)
                        NoiseStart = NoiseStart - datetime.timedelta(hours = NoiseStart.hour, minutes = NoiseStart.minute, seconds=NoiseStart.second) #rounds down to the nearest day
                        NoiseEnd = UTCDateTime.strptime(ev,dateformat)
                        NoiseEnd = NoiseEnd - datetime.timedelta(hours = NoiseEnd.hour, minutes = NoiseEnd.minute, seconds=NoiseEnd.second) #rounds down to the nearest day
                        args = [staquery_output,'--start={}'.format(NoiseStart), '--end={}'.format(NoiseEnd)]
                        # sys.stdout.flush()
                        atacr_download_data.main(atacr_download_data.get_daylong_arguments(args))
        print(' ')
        print('----Noise Download Complete----')
        sys.stdout = original

def DailySpectra(catalog,netsta_names=None,extra_flags = '-O --figQC --figAverage --figCoh --save-fig',subfolder=None):
        # sys.stdout.flush()
        datafolder = './Data/'
        staquery_output = datafolder + 'sta_query.pkl'
        if subfolder is not None:
                logoutput = subfolder + '_Step_4_7_QCSpectra_logfile.log'
        else:
                logoutput = '_Step_4_7_QCSpectra_logfile.log'
        log_fout = datafolder + logoutput
        SpecStart = catalog.Start.min().strftime("%Y-%m-%d, %H:%M:%S")
        SpecEnd = catalog.End.max().strftime("%Y-%m-%d, %H:%M:%S")
        args = [staquery_output]
        [args.append(flg) for flg in extra_flags.split()]
        [args.append(flg) for flg in ['--start={}'.format(SpecStart),'--end={}'.format(SpecEnd)]]
        # with open(log_fout, 'w') as sys.stdout:
        original = sys.stdout
        sys.stdout = open(log_fout,'w+')
        print('----Begin Daily Spectra----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.io.build_staquery(catalog,staquery_output)
        atacr_daily_spectra.main(atacr_daily_spectra.get_dailyspec_arguments(args))
        print(' ')
        print('----Daily Spectra Complete----')
        sys.stdout = original

def CleanSpectra(catalog,netsta_names=None,extra_flags = '-O --figQC --figAverage --figCoh --figCross --save-fig',subfolder=None):
        # sys.stdout.flush()
        datafolder = './Data/'
        staquery_output = datafolder + 'sta_query.pkl'
        if subfolder is not None:
                logoutput = subfolder + '_Step_5_7_CleanSpectra_logfile.log'
        else:
                logoutput = '_Step_5_7_CleanSpectra_logfile.log'
        log_fout = datafolder + logoutput
        SpecStart = catalog.Start.min().strftime("%Y-%m-%d, %H:%M:%S")
        SpecEnd = catalog.End.max().strftime("%Y-%m-%d, %H:%M:%S")
        args = [staquery_output]
        [args.append(flg) for flg in extra_flags.split()]
        [args.append(flg) for flg in ['--start={}'.format(SpecStart),'--end={}'.format(SpecEnd)]]
        # with open(log_fout, 'w') as sys.stdout:
        original = sys.stdout
        sys.stdout = open(log_fout,'w+')
        print('----Begin Clean Spectra----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.io.build_staquery(catalog,staquery_output)
        atacr_clean_spectra.main(atacr_clean_spectra.get_cleanspec_arguments(args))
        print(' ')
        print('----Clean Spectra Complete----')
        sys.stdout = original

def TransferFunctions(catalog,netsta_names=None,extra_flags = '-O --figTF --save-fig',subfolder=None):
        # sys.stdout.flush()
        datafolder = './Data/'
        staquery_output = datafolder + 'sta_query.pkl'
        if subfolder is not None:
                logoutput = subfolder + '_Step_6_7_CalcTFs_logfile.log'
        else:
                logoutput = '_Step_6_7_CalcTFs_logfile.log'
        log_fout = datafolder + logoutput
        args = [staquery_output]
        [args.append(flg) for flg in extra_flags.split()]
        # with open(log_fout, 'w') as sys.stdout:
        original = sys.stdout
        sys.stdout = open(log_fout,'w+')
        print('----Begin Building Transfer Functions----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.io.build_staquery(catalog,staquery_output)
        atacr_transfer_functions.main(atacr_transfer_functions.get_transfer_arguments(args))

        print(' ')
        print('----Building Transfer Functions Complete----')
        sys.stdout = original

def CorrectEvents(catalog,netsta_names=None,extra_flags = '--figRaw --figClean --save-fig',subfolder=None):
        # sys.stdout.flush()
        datafolder = './Data/'
        staquery_output = datafolder + 'sta_query.pkl'
        if subfolder is not None:
                logoutput = subfolder + '_Step_7_7_CorrectEvents_logfile.log'
        else:
                logoutput = '_Step_7_7_CorrectEvents_logfile.log'
        log_fout = datafolder + logoutput
        args = [staquery_output]
        [args.append(flg) for flg in extra_flags.split()]
        # with open(log_fout, 'w') as sys.stdout:
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.io.build_staquery(catalog,staquery_output)
        original = sys.stdout
        sys.stdout = open(log_fout,'w+')
        print('----Begin Correcting Eventss----')
        atacr_correct_event.main(atacr_correct_event.get_correct_arguments(args))
        print(' ')
        print('----Correcting Events Complete----')
        sys.stdout = original



def Run_ATaCR(catalog, STEPS=[2,3,4,5,6,7], netsta_names=None, Minmag=6.3, Maxmag=6.7, limit=1000, pre_event_min_aperture=1, pre_event_day_aperture=30, dailyspectra_flags='-O --figQC --figAverage --figCoh --save-fig', cleanspectra_flags='-O --figQC --figAverage --figCoh --figCross --save-fig', tf_flags='-O --figTF --save-fig', correctevents_flags='--figRaw --figClean --save-fig',subfolder=None):

        # STEPS = [1,2,3,4,5,6,7] #Absolutely every step - Downloading adds an hour or more to the process
        # STEPS = [2,3] #Everything but the download steps - About 4min for six stations.
        # STEPS = [4,5,6,7] #Everything but the download steps - About 4min for six stations.

        if 1 in STEPS:
                print('Step 1/7 - BEGIN: Station Metadata')
                # C='?H?' #channels
                # !query_fdsn_stdb -N {','.join(N)} -C '{C}' -S {','.join(S)} ./Data/sta_query> ./Data/Step_1_7_StationMeta_logfile.log
                ObsQA.io.build_staquery(catalog,staquery_output = './Data/sta_query.pkl',subfolder=subfolder)
                print('Step 1/7 - COMPLETE: Station Metadata')
        if 2 in STEPS:
                print('Step 2/7 - BEGIN: Download Event Data')
                ObsQA.io.DownloadEvents(catalog,netsta_names=netsta_names,Minmag=Minmag,Maxmag=Maxmag,limit=limit,pre_event_min_aperture=pre_event_min_aperture,subfolder=subfolder)
                print('Step 2/7 - COMPLETE: Download Event Data')
        if 3 in STEPS:
                print('Step 3/7 - BEGIN: Download Day Data')
                ObsQA.io.DownloadNoise(catalog,netsta_names=netsta_names,pre_event_day_aperture=pre_event_day_aperture,subfolder=subfolder)
                print('Step 3/7 - COMPLETE: Download Day Data')

        if 4 in STEPS:
                print('Step 4/7 - BEGIN: Quality Control Noise Data')
                ObsQA.io.DailySpectra(catalog,netsta_names=netsta_names,extra_flags=dailyspectra_flags,subfolder=subfolder)
                print('Step 4/7 - COMPLETE: Quality Control Noise Data')
        if 5 in STEPS:
                print('Step 5/7 - BEGIN: Spectral Average of Noise Data')
                # !atacr_clean_spectra -O --figQC --figAverage --figCoh --figCross --save-fig --start='{SpecStart}' --end='{SpecEnd}' ./Data/sta_query.pkl> ./Data/Step_5_7_CleanSpectra_logfile.log
                ObsQA.io.CleanSpectra(catalog,netsta_names=netsta_names,extra_flags=cleanspectra_flags,subfolder=subfolder)
                print('Step 5/7 - COMPLETE: Spectral Average of Noise Data')
        if 6 in STEPS:
                print('Step 6/7 - BEGIN: Calculate Transfer Functions')
                # !atacr_transfer_functions -O --figTF --save-fig ./Data/sta_query.pkl> ./Data/Step_6_7_CalcTFs_logfile.log
                ObsQA.io.TransferFunctions(catalog,netsta_names=netsta_names,extra_flags=tf_flags,subfolder=subfolder)
                print('Step 6/7 - COMPLETE: Calculate Transfer Functions')
        if 7 in STEPS:
                print('Step 7/7 - BEGIN: Correct Event Data')
                # !atacr_correct_event --figRaw --figClean --save-fig ./Data/sta_query.pkl> ./Data/Step_7_7_CorrectEvents_logfile.log
                ObsQA.io.CorrectEvents(catalog,netsta_names=netsta_names,extra_flags=correctevents_flags,subfolder=subfolder)
                print('Step 7/7 - COMPLETE: Correct Event Data')

# eof