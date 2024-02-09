from ObsQA import *
import ObsQA as ob
from ObsQA.classes import OBSMetrics
from scipy.stats import norm
from cmcrameri import cm
import matplotlib.colors as mcolors

def OBS_Generator(catalog,datafolder,events_folder = 'EVENTS',avg='STA',tf='ZP-21',evireq=None,sta=None, net=None,hps_win_length=200,width=None):
        # datafolder = ATaCR_Py_DataFolder['Py_CorrectedTraces']
        nsta = len(catalog)
        cat = catalog.copy()
        staind = [i for i in range(len(cat))]
        if sta is not None:
            fi = np.atleast_1d(np.in1d(cat.Station,np.array(sta).squeeze()))
            cat = cat[fi].copy()
            staind = staind[np.squeeze(np.where(fi))]
        if net is not None:
            fi = np.in1d(cat.Network,np.array(net))
            cat = cat[fi].copy()
            staind = staind[np.squeeze(np.where(fi))]
        staind = np.atleast_1d(np.array(staind))
        for stai,Station in zip(staind,cat.iloc):
                StaName = Station.StaName
                Metrics = dict()
                evts = list(Station.Events)
                evind = [i for i in range(len(evts))]
                if evireq is not None:
                    evts = evts[evireq]
                    evind = evind[evireq]
                nevents = len(evts)
                # folder = datafolder + '/' + StaName + '/Corrected'
                try:
                    Noise = get_Noise(datafolder,Station.Network,Station.Station,avg='STA',update=True)
                    Metrics['Noise'] = ObsQA.classes.OBSMetrics(csd=Noise.Noise[0],f=Noise.Noise[0].f)
                    Metrics['StaNoise'] = Noise
                except:
                    print('Noise files not found')
                for evi,Event in zip(evind,evts):
                        Event = Event + '  ' + 'm' + str(Station.Magnitude_mw[evi]) + '  ' + str(Station.Depth_KM[evi]) + 'km'
                        Notify = StaName + ' | ' + str(stai+1) + '/' + str(nsta) + ' | ' + str(evi+1) + '/' + str(nevents) + ' : ' + Station.Events[evi]
                        Comps = dict()
                        try:
                                d = ob.io.Get_ATaCR_CorrectedEvents(datafolder + '/' + events_folder,[Station.Events[evi]],Station.Network,Station.Station,tfavg=avg.lower())
                                Comps['d'] = d
                        except:
                                # print(Notify + ' : ' + 'No File')
                                continue
                        # MetricKeys = list(d.Corrected[0].keys())
                        Comps['RawP'],Comps['RawZ'],Comps['PostATACR'],Comps['PostHPS'],Comps['PostBoth'],Comps['PostHPS_H1'],Comps['PostHPS_H2'] = get_ATACR_HPS_Comp(d,atacr_TFs_used=tf, win_length=hps_win_length,width=width)
                        Metrics['Raw'] = OBSMetrics(         tr1=d.Raw[0]['tr1'].copy(),       tr2=d.Raw[0]['tr2'].copy(),       trZ=d.Raw[0]['trZ'].copy(),       trP=d.Raw[0]['trP'].copy())
                        Metrics['ATaCR'] = OBSMetrics(       tr1=d.Raw[0]['tr1'].copy(),       tr2=d.Raw[0]['tr2'].copy(),       trZ=d.Corrected[0][tf].copy(),    trP=d.Raw[0]['trP'].copy()) / Metrics['Raw'].copy()
                        Metrics['HPS'] = OBSMetrics(         tr1=Comps['PostHPS_H1'].copy(),   tr2=Comps['PostHPS_H2'].copy(),   trZ=Comps['PostHPS'].copy(),      trP=d.Raw[0]['trP'].copy()) / Metrics['Raw'].copy()
                        Metrics['ATaCR_HPS'] = OBSMetrics(   tr1=Comps['PostHPS_H1'].copy(),   tr2=Comps['PostHPS_H2'].copy(),   trZ=Comps['PostBoth'].copy(),     trP=d.Raw[0]['trP'].copy()) / Metrics['Raw'].copy()
                        # Metrics['Raw w/ Noise'] = Metrics['Raw'].copy() / Metrics['Noise'].copy()
                        Metrics['ATaCR w/ HPS'] = Metrics['ATaCR'].copy() / Metrics['HPS'].copy()
                        Metrics['ATaCR_HPS w/ HPS'] = Metrics['ATaCR_HPS'].copy() / Metrics['HPS'].copy()
                        Metrics['ATaCR_HPS w/ ATaCR'] = Metrics['ATaCR_HPS'].copy() / Metrics['ATaCR'].copy()
                        print(Notify)
                        yield (Event,Station,Metrics,Comps)

def get_ATACR_HPS_Comp(d,atacr_TFs_used='ZP-21', win_length=200,width=None):
    RawP = d.Raw[0]['trP'].copy()
    RawZ = d.Raw[0]['trZ'].copy()
    PostATACR = d.Corrected[0][atacr_TFs_used].copy()
    PostATACR.stats.location = 'ATaCR (' + PostATACR.stats.location + ')'
    PostHPS, spectrograms = noisecut(RawZ.copy(), ret_spectrograms=True,win_length=win_length,width=width)
    PostBoth, spectrograms = noisecut(PostATACR.copy(), ret_spectrograms=True,win_length=win_length,width=width)
    PostHPS_H1,spectrograms = noisecut(d.Raw[0]['tr1'].copy(), ret_spectrograms=True,win_length=win_length,width=width)
    PostHPS_H2,spectrograms = noisecut(d.Raw[0]['tr2'].copy(), ret_spectrograms=True,win_length=win_length,width=width)
    # PostHPS_P,spectrograms = noisecut(d.Raw[0]['trP'].copy(), ret_spectrograms=True,win_length=win_length,width=width)
    return RawP,RawZ,PostATACR,PostHPS,PostBoth,PostHPS_H1,PostHPS_H2

def get_obs_cmaps(light_min=0.4,light_max=0.9,extrema_n=4):
    cmaps = dict()
    raw_palette = 'grayC_categorical'
    atacr_palette = 'bilbao_categorical'
    nc_palette = 'devon_categorical'
    both_palette = 'acton_categorical'
    comp_palette = 'bamako_categorical'
    for k,p in zip(['raw','atacr','nc','both','comp'],[raw_palette,atacr_palette,nc_palette,both_palette,comp_palette]):
            colors = cm.cmaps[p].colors
            colors = cm.cmaps[p].colors[((cm.cmaps[p].colors**2).sum(axis=1)**(0.5))<light_max]
            colors = colors[(colors**2).sum(axis=1)**(0.5)>light_min]
            if k=='atacr':
                    colors[0][0]  = colors[0][0]**(0.5**(extrema_n+1))
                    for i,c in enumerate(colors):
                        c[0] = c[0]**(0.5**extrema_n)
                        colors[i] = c
            if k=='nc':
                    colors[0][2]  = colors[0][2]**(0.5**(extrema_n+1))
                    for i,c in enumerate(colors):
                        c[2] = c[2]**(0.5**extrema_n)
                        colors[i] = c
            if k=='raw':
                    # pass
                    colors[0] = [0,0,0]
                    # colors[0] = colors[0]**(2**extrema_n)
                    # for i,c in enumerate(colors):
                    #     c = c**(2**extrema_n)
                    #     colors[i] = c
            if k=='both':
                    colors[0][0]  = colors[0][0]**(0.5**(extrema_n+1))
                    colors[0][2]  = colors[0][2]**(0.5**(extrema_n+1))
                    for i,c in enumerate(colors):
                        c[0] = c[0]**(0.5**extrema_n)
                        c[2] = c[2]**(0.5**extrema_n)
                        colors[i] = c
            if k=='comp':
                    colors[0][1]  = colors[0][1]**(0.5**(extrema_n+1))
                    for i,c in enumerate(colors):
                        c[1] = c[1]**(0.5**extrema_n)
                        colors[i] = c
            colors = [mcolors.to_hex(c) for c in colors]
            cmaps[k] = colors
    return cmaps

def log_smoothing(x,y):
    # Peak Log-Smoothed trendline.
    # Math is inefficient but it works. Needs to be replaced with np.ceil(np.log10(Nyq/NFFT))
    dw = x[np.log10(x)<=np.ceil(np.log10(x[1]))].size
    dx = x[np.log10(x)<=np.ceil(np.log10(x[1]))].max()
    kernel_size = int(np.round((dx/x[-1])*x.size))
    kernel = np.ones(kernel_size) / kernel_size
    y = np.convolve(y, kernel, mode='valid')
    x = x[np.arange(0,y.size,dw)]
    y = y[np.arange(0,y.size,dw)]
    return x,y


def save_tight(filename,fig=None,dpi=200):
    # Saves figure to PDF with no margins. Do not modify
    # plt.gca().set_axis_off()
    plt.subplots_adjust(top = 1, bottom = 0.0, right = 1, left = 0,hspace = 0.07, wspace = 0.03)
    plt.margins(0.1,0.1)
    # plt.gca().xaxis.set_major_locator(plt.NullLocator())
    # plt.gca().yaxis.set_major_locator(plt.NullLocator())
    if fig is None:
        plt.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi)
    else:
        fig.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi)
    return 'Complete'