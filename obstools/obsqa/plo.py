def ML_fig_event_corrected(PreProcEvent,CorrectedEvent,evstream=None, TF_list=None, fmin=1./150., fmax=2.,prefix = '',yes_filter=True,ylon=True,scale='linear',yhard=None):
        """
        Adapted from Python ATaCR code. Refactored for efficiency with a few more options. -CH-8/17/23
        Function to plot the corrected vertical component seismograms.
        Parameters
        ----------
        evstream : :class:`~obtsools.classes.EventStream`
                Container for the event stream data
        Tf_list : list
                List of Dictionary elements of transfer functions used
                for plotting the corrected vertical component.
        """
        keys = ['Z1','Z2-1','ZP-21','ZH','ZP-H','ZP']
        # Unpack vertical trace and filter
        preproc_itr1 = np.where((PreProcEvent['channel'].squeeze().str.find('1')>0))[0][0]
        preproc_itr2 = np.where((PreProcEvent['channel'].squeeze().str.find('2')>0))[0][0]
        preproc_itrZ = np.where((PreProcEvent['channel'].squeeze().str.find('Z')>0))[0][0]
        preproc_itrP = np.where(~((PreProcEvent['channel'].squeeze().str.find('1')>0) + (PreProcEvent['channel'].squeeze().str.find('2')>0) + (PreProcEvent['channel'].squeeze().str.find('Z')>0)))[0][0]
        trZ = PreProcEvent.iloc[preproc_itrZ]['Trace']
        trZ.filter(
                'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
        sr = trZ.stats.sampling_rate
        taxis = np.arange(0., trZ.stats.npts/sr, 1./sr)
        if ylon:
                yl = [-1e-4,1e-4]
        if yhard is not None:
                yl = yhard
        plt.figure(figsize=(8, 8))
        eventtime = CorrectedEvent.eventid.squeeze().to_list()[0]
        evtstamp = str(UTCDateTime(eventtime).year) + '.' + str(UTCDateTime(eventtime).julday) + '.' + str(UTCDateTime(eventtime).hour) + '.' + str(UTCDateTime(eventtime).minute)
        i = 611
        for key in keys:
                plt.subplot(i)
                if scale=='symlog':
                        plt.plot(taxis, np.abs(trZ.data), 'lightgray', lw=0.5)
                else:
                        plt.plot(taxis, trZ.data, 'lightgray', lw=0.5)
                plt.yscale(scale)
                if TF_list[key]:
                        Corrtrace = CorrectedEvent.iloc[np.where(~CorrectedEvent.label.squeeze().str.find(key))[0][0]]
                        tr = Trace(data=Corrtrace['Trace'].data,header=Corrtrace['Trace'].stats)
                        if yes_filter:
                                tr = tr.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
                        if scale=='symlog':
                                plt.plot(taxis, np.abs(tr.data), 'k', lw=0.5)
                                plt.plot(taxis, np.abs(tr.data), 'k', lw=0.5)
                        else:
                                plt.plot(taxis, tr.data, 'k', lw=0.5)
                                plt.plot(taxis, tr.data, 'k', lw=0.5)
                        plt.yscale(scale)
                plt.title(prefix + Corrtrace.network + '.' + Corrtrace.station + ' ' + evtstamp + ': ' + key, fontdict={'fontsize': 8})
                if scale=='linear':
                        plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,scilimits=(-3, 3))
                plt.xlim((0., trZ.stats.npts/sr))
                if ylon:
                        plt.ylim(yl)
                if yhard is not None:
                        plt.ylim(yl)
                i+=1
        plt.xlabel('Time since earthquake (sec)')
        plt.tight_layout()
        return plt

def Py_fig_event_corrected(evstream, TF_list, fmin=1./150., fmax=2.,prefix=''):
        """
        Adapted ATaCR plot code with a few more options for versatility - CH-8/15/23
        Function to plot the corrected vertical component seismograms.
        Parameters
        ----------
        evstream : :class:`~obtsools.classes.EventStream`
                Container for the event stream data
        Tf_list : list
                List of Dictionary elements of transfer functions used
                for plotting the corrected vertical component.
        """

        # Unpack vertical trace and filter
        trZ = evstream.trZ.copy()
        trZ.filter(
                'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
        sr = trZ.stats.sampling_rate
        taxis = np.arange(0., trZ.stats.npts/sr, 1./sr)
        plt.figure(figsize=(8, 8))
        plt.subplot(611)
        plt.plot(
                taxis, trZ.data, 'lightgray', lw=0.5)
        if TF_list['Z1']:
                tr = Trace(
                data=evstream.correct['Z1'],
                header=trZ.stats).filter(
                'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
                plt.plot(taxis, tr.data, 'k', lw=0.5)
        plt.title(prefix + evstream.key + ' ' + evstream.tstamp +
                ': Z1', fontdict={'fontsize': 8})
        plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                                scilimits=(-3, 3))
        plt.xlim((0., trZ.stats.npts/sr))
        plt.subplot(612)
        plt.plot(taxis, trZ.data, 'lightgray', lw=0.5)
        if TF_list['Z2-1']:
                tr = Trace(
                data=evstream.correct['Z2-1'],
                header=trZ.stats).filter(
                'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
                plt.plot(taxis, tr.data, 'k', lw=0.5)
        plt.title(prefix + evstream.tstamp + ': Z2-1', fontdict={'fontsize': 8})
        plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                                scilimits=(-3, 3))
        plt.xlim((0., trZ.stats.npts/sr))
        plt.subplot(613)
        plt.plot(
                taxis, trZ.data, 'lightgray', lw=0.5)
        if TF_list['ZP-21']:
                tr = Trace(
                data=evstream.correct['ZP-21'],
                header=trZ.stats).filter(
                'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
                plt.plot(taxis, tr.data, 'k', lw=0.5)
        plt.title(prefix + evstream.tstamp + ': ZP-21', fontdict={'fontsize': 8})
        plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                                scilimits=(-3, 3))
        plt.xlim((0., trZ.stats.npts/sr))
        plt.subplot(614)
        plt.plot(
                taxis, trZ.data, 'lightgray', lw=0.5)
        if TF_list['ZH']:
                tr = Trace(
                data=evstream.correct['ZH'],
                header=trZ.stats).filter(
                'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
                plt.plot(taxis, tr.data, 'k', lw=0.5)
        plt.title(prefix + evstream.tstamp + ': ZH', fontdict={'fontsize': 8})
        plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                                scilimits=(-3, 3))
        plt.xlim((0., trZ.stats.npts/sr))
        plt.subplot(615)
        plt.plot(taxis, trZ.data, 'lightgray', lw=0.5)
        if TF_list['ZP-H']:
                tr = Trace(
                data=evstream.correct['ZP-H'],
                header=trZ.stats).filter(
                'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
                plt.plot(taxis, tr.data, 'k', lw=0.5)
        plt.title(prefix + evstream.tstamp + ': ZP-H', fontdict={'fontsize': 8})
        plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                                scilimits=(-3, 3))
        plt.xlim((0., trZ.stats.npts/sr))

        plt.subplot(616)
        plt.plot(
                taxis, trZ.data, 'lightgray', lw=0.5)
        if TF_list['ZP']:
                tr = Trace(
                data=evstream.correct['ZP'],
                header=trZ.stats).filter(
                'bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
                plt.plot(taxis, tr.data, 'k', lw=0.5)
        plt.title(prefix + evstream.tstamp + ': ZP', fontdict={'fontsize': 8})
        plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,
                                scilimits=(-3, 3))
        plt.xlim((0., trZ.stats.npts/sr))
        plt.xlabel('Time since earthquake (sec)')
        plt.tight_layout()
        return plt

def PlotPh(input):
        '''
        Takes a single obstools.atacr.classes.StaNoise object and plots spectral PHASE DISPARITY between all component pairs contained within.
        '''
        ph_map = {0:'c1Z',1:'c2Z',2:'cHZ',3:'cZP',4:'c12',5:'c1P',6:'c2P'} #DO NOT MODIFY!
        if not isinstance(input,list):
                input = [input]
        cols=len(input)
        plt.figure(figsize=( (16+16*(0.25*(cols-1))),2*len(ph_map) ))
        for M, io in zip(input,range(cols)):
                for i in range(len(ph_map)):
                        a = ph_map[i]
                        ph = Phase(M.__dict__[a])
                        ttl = M.key + ' - Spectral phase disparity' + ' - ' + str(np.sum(M.nwins)) + ' windows over ' + str(str(ph.shape[1])) + ' days'
                        ttlb = '\nbetween ' + a[1] + ' and ' + a[2]
                        for e in range(ph.shape[1]):
                                plt.subplot(len(ph_map),cols,cols*(i+io)+1-io)
                                plt.subplots_adjust(bottom=0.1,left=0.05 ,right=0.7, top=0.9)
                                plt.scatter(M.f,ph[:,e],c='k',s=0.5)
                                plt.xscale('log')
                                if io==0:
                                        plt.ylabel('Phase (Degrees)')
                                if i==0:
                                        plt.title(ttl + ttlb)
                                else:
                                        plt.title(ttlb)
                                if i==(len(ad_map)-1):
                                        plt.xlabel('Frequency (Hz)')
                                if i<(len(ad_map)-1):
                                        plt.gca().set_xticklabels([])

def PlotAd(input):
        '''
        Takes a single obstools.atacr.classes.StaNoise object and plots spectral ADMITTANCE between all component pairs contained within.
        '''
        ad_map = {0:['c1Z','c11'],1:['c2Z','c22'],2:['cHZ','cHH'],3:['cZP','cPP'],4:['c12','c11'],5:['c1P','c11'],6:['c2P','c22']} #DO NOT MODIFY!
        if not isinstance(input,list):
                input = [input]
        cols=len(input)
        plt.figure(figsize=( (16+16*(0.25*(cols-1))),2*len(ad_map) ))
        for M, io in zip(input,range(cols)):
                for i in range(len(ad_map)):
                        a,b = ad_map[i]
                        ad = Admittance(M.__dict__[a],M.__dict__[b])
                        ttl = M.key + ' - Spectral admittance' + ' - ' + str(np.sum(M.nwins)) + ' windows over ' + str(str(ad.shape[1])) + ' days'
                        ttlb = '\nbetween '+ a[1] + ' and ' + a[2]
                        for e in range(ad.shape[1]):
                                plt.subplot(len(ad_map),cols,cols*(i+io)+1-io)
                                plt.subplots_adjust(bottom=0.1,left=0.05 ,right=0.7, top=0.9)
                                plt.plot(np.atleast_2d(M.f).T,ad[:,e].T,c='k',linewidth=0.1)
                                plt.xscale('log')
                                plt.yscale('log')
                                if io==0:
                                        plt.ylabel('Admittance')
                                if i==0:
                                        plt.title(ttl + ttlb)
                                else:
                                        plt.title(ttlb)
                                if i==(len(ad_map)-1):
                                        plt.xlabel('Frequency (Hz)')
                                if i<(len(ad_map)-1):
                                        plt.gca().set_xticklabels([])

def PlotCOH(input):
        '''
        Takes a single obstools.atacr.classes.StaNoise object and plots spectral COHERENCE between all component pairs contained within.
        '''
        coh_map = {0:['c1Z','c11','cZZ'],1:['c2Z','c22','cZZ'],2:['cHZ','cHH','cZZ'],3:['cZP','cPP','cZZ'],4:['c12','c11','c22'],5:['c1P','c11','cPP'],6:['c2P','c22','cPP']} #DO NOT MODIFY!
        if not isinstance(input,list):
                input = [input]
        cols=len(input)
        plt.figure(figsize=( (16+16*(0.25*(cols-1))),2*len(coh_map) ))
        for M, io in zip(input,range(cols)):
                for i in range(len(coh_map)):
                        a,b,c = coh_map[i]
                        coh = Coherence(M.__dict__[a],M.__dict__[b],M.__dict__[c])
                        ttl = M.key + ' - Spectral coherence' + ' - ' + str(np.sum(M.nwins)) + ' windows over ' + str(str(coh.shape[1])) + ' days'
                        ttlb = '\nbetween '+ a[1] + ' and ' + a[2]
                        for e in range(coh.shape[1]):
                                plt.subplot(len(coh_map),cols,cols*(i+io)+1-io)
                                plt.subplots_adjust(bottom=0.1,left=0.05 ,right=0.7, top=0.9)
                                plt.plot(np.atleast_2d(M.f).T,coh[:,e],c='k',linewidth=0.1)
                                plt.xscale('log')
                                if io==0:
                                        plt.ylabel('Coherence')
                                if i==0:
                                        plt.title(ttl + ttlb)
                                else:
                                        plt.title(ttlb)
                                if i==(len(coh_map)-1):
                                        plt.xlabel('Frequency (Hz)')
                                if i<(len(coh_map)-1):
                                        plt.gca().set_xticklabels([])
# eof