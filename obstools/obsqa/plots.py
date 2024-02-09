from ObsQA.imports import *
#---I don't have a spot for these in the package yet. Still workging on the PSD testing::::
# import scipy.stats as stats
# import scipy as sc
# from numpy.lib.stride_tricks import sliding_window_view
#---Not used in my QA. Just have it here for reference:::
# from obstools.atacr.plotting import fig_QC, fig_average, fig_av_cross, fig_coh_ph, fig_TF, fig_comply, fig_event_raw, fig_event_corrected

def CompPlot(x,pre,post,group_title=['',0],style='plot',title = '',xlabel='',ylabel='',xscale='linear',yscale='linear',alpha=(0.2,1.0),color=('k','k'),figsize=(10,10),share='row',s=[0.5,0.2],bg=None,bg_label = None,pre_label=None,post_label=None,lgloc='upper left'):
        #  x: an ndarray, array, list,  or list of lists each containing the x-array used in each row of the subplots
        #  pre and post: An ndarray, array,list, or list of lists containing the pre/post y-values to be plotted
        #  xlabel and ylabel : Either a str or list objects such that the ylabel is the same
        #  length as the x list and likewise for xlabel with the pre/post list lengths.
        #  xscale and yscale: a str or list of str of the same lengths of columns and rows, respectively
        #  style: 'scatter' or 'plot'
        #  alpha[0]/alpha[1]: float
        #  figsize: tuple
        pre = np.atleast_2d(np.array(pre)).squeeze()
        post = np.atleast_2d(np.array(post)).squeeze()
        if len(pre.shape)<3:
                pre = np.atleast_3d(np.array(pre))
                post = np.atleast_3d(np.array(post))
                pre = np.moveaxis(pre,[0,1,2],[1,2,0])
                post = np.moveaxis(post,[0,1,2],[1,2,0])
        x = np.atleast_2d(np.array(x))
        if x.shape[0] is not pre.shape[0]:
                x = x.repeat(3,axis=0)
        if isinstance(title,list):
                title = np.atleast_2d(np.array(title))
                if title.size==pre[:,:,0].size:
                        title = title.reshape(pre[:,:,0].shape)
        if not isinstance(xscale,list):
                xscale = np.array([xscale for  _ in range(x.shape[0])]).reshape(-1)
                yscale = np.array([yscale for  _ in range(x.shape[0])]).reshape(-1)
        if not isinstance(xlabel,list):
                xlabel = np.array([xlabel for  _ in range(x.shape[0])]).reshape(-1)
        if not isinstance(ylabel,list):
                ylabel = np.array([ylabel for  _ in range(x.shape[0])]).reshape(-1)
        if not isinstance(style,list):
                style = np.array([style for  _ in range(x.shape[0])]).reshape(-1)
        if not isinstance(s,list):
                s = np.array([s for  _ in range(2)]).reshape(-1)
        fig,ax = plt.subplots(x.shape[0],pre.shape[len(pre.shape)-2],figsize=figsize,sharex=share,sharey=share)
        ax = np.atleast_2d(ax)
        fig.tight_layout()
        for xi,x_plot in enumerate(x):
                x_plot = x_plot.reshape(-1)
                pre_row = np.atleast_2d(pre[xi,:])
                post_row = np.atleast_2d(post[xi,:])
                axlims = [pre_row.min(),pre_row.max(),post_row.min(),post_row.max()]
                axlims = np.round([np.min(axlims),np.max(axlims)])
                for coli,(current_axis,y_pre,y_post) in enumerate(zip(ax[xi],pre_row,post_row)):
                        y_pre = y_pre.reshape(-1)
                        y_post = y_post.reshape(-1)
                        if bg is not None:
                                bg_path = current_axis.scatter(bg[xi][0],bg[xi][1],alpha=0.04,s=np.min(s)/2,edgecolors=None,c='grey')
                                if bg_label is not None:
                                        bg_path.set_label(bg_label)
                        if style[xi].lower()=='scatter':
                                prepath = current_axis.scatter(x_plot,y_pre,alpha=alpha[0],c=color[0],s=s[0]/1,edgecolors=None)
                                postpath = current_axis.scatter(x_plot,y_post,alpha=alpha[1],c=color[1],s=s[1]/1000,edgecolors=None)
                        elif style[xi].lower()=='plot':
                                prepath = current_axis.plot(x_plot,y_pre,alpha=alpha[0],c=color[0],linewidth=s[0])
                                postpath = current_axis.plot(x_plot,y_post,alpha=alpha[1],c=color[1],linewidth=s[1])
                        if (xi==0) & (coli==0):
                                if pre_label is not None:
                                        prepath.set_label(pre_label)
                                        lg = current_axis.legend(loc=lgloc,fontsize=9)
                                if post_label is not None:
                                        postpath.set_label(post_label)
                                        lg = current_axis.legend(loc=lgloc,fontsize=9)
                                for i in range(len(lg.legendHandles)):
                                        lg.legendHandles[i]._sizes = np.array([20])
                                        lg.legendHandles[i]._alpha = 1.0
                                        # lg.legendHandles[i]._sizes = np.array([100])
                        # current_axis.autoscale(enable=True, axis='both', tight=True)
                        # current_axis.set_xlim(x_plot.min(),x_plot.max())
                        # current_axis.set_ylim(axlims[0],axlims[1])
                        current_axis.set_xscale(xscale[xi])
                        current_axis.set_yscale(yscale[xi])
                        if (xi==0) & (coli==1) & isinstance(title,str):
                                current_axis.set_title(title)
                        elif isinstance(title, (np.ndarray, np.generic)):
                                if (title.shape[0]-1)>=xi:
                                        if (xi==0) & (coli==group_title[1]):
                                                current_axis.set_title(group_title[0] + '\n' + title[xi,coli],fontweight='bold')
                                        else:
                                                current_axis.set_title(title[xi,coli],fontweight='bold')
                                elif xi==0:
                                        current_axis.set_title(title[0,coli],fontweight='bold')
                        if coli==0:
                                current_axis.set_ylabel(ylabel[xi],fontweight='bold',fontsize=12)
                        if xi==(x.shape[0]-1):
                                current_axis.set_xlabel(xlabel[0])
                        else:
                                current_axis.set_xticklabels('')
        plt.tight_layout()

def plotPrePost(Pre,Post,ysig=2,fmin=None,fmax=None,yl=None,xl=None,title=None,ev_llaz=None):
        if not isinstance(Pre,list):
                Pre = [Pre]
        if not isinstance(Post,list):
                Post = [Post]
        alpha =[0.2,1.0]
        nrows = len(Pre)
        plt.figure(figsize=(20,3*nrows))
        for i in range(nrows):
                plt.subplot(nrows,1,i+1)
                PreA = Pre[i]
                PostB = Post[i]
                ts = [PreA.times(),PostB.times()]
                trs = [PreA,PostB]
                for a in range(len(alpha)):
                        t = ts[a]
                        y = trs[a].copy()
                        if fmin is None and fmax is not None:
                                y.filter('lowpass', freq=fmax,corners=4,zerophase=True)
                        elif fmin is not None and fmax is None:
                                y.filter('highpass', freq=fmin, corners=4, zerophase=True)
                        elif fmin is not None and fmax is not None:
                                y.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)
                        plt.plot(t,y.data,alpha=alpha[a],color='k')
                        if ysig is not None:
                                plt.ylim(np.mean(y.data)-ysig*np.std(y.data),np.mean(y.data)+ysig*np.std(y.data))
                        elif yl is not None:
                                plt.ylim([-yl,yl])
                        else:
                                if a==1:
                                        plt.ylim([np.min(y.data),np.max(y.data)])
                        if xl is not None:
                                plt.xlim(xl)
                        else:
                                plt.xlim(t[0],t[-1])
                        if a==1:
                                if ev_llaz is not None:
                                        sta_llaz = [y.stats.sac.stla,y.stats.sac.stlo,np.abs(y.stats.sac.stel)]
                                        arrivals = ObsQA.io.get_arrivals(sta_llaz,ev_llaz)
                                        ObsQA.plots.overplot_taup(arrivals)
                if (i+1)==nrows:
                        plt.xlabel('Seconds')
                else:
                        plt.gca().set_xticklabels('')
                # if i==0:
                if i==0 and title is not None:
                        plt.title(title + '\n' + Pre[i].stats.location + ' >> ' + Post[i].stats.location)
                else:
                        plt.title(Pre[i].stats.location + ' >> ' + Post[i].stats.location)
                plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))


def overplot_taup(arrivals):
        ii=4
        ang = 30
        for a in arrivals:
                ii+=1
                y = plt.gca().get_ylim()
                xl = plt.gca().get_xlim()
                x = [a[1],a[1]]
                alg = ['enter','left','right']
                e = -1
                if a[0][-1].upper()=='S':
                        c = 'b'
                else:
                        c = 'r'
                if (len([m.start() for m in re.finditer('S',a[0].upper())])<=1) or (len([m.start() for m in re.finditer('I',a[0].upper())])>=1) or (len([m.start() for m in re.finditer('C',a[0].upper())])<1):
                        if x[0]<np.max(xl):
                                if len(a[0])<2:
                                        plt.plot(x,y,c,linewidth=1,linestyle='-.')
                                        if a[0][-1].upper()=='S':
                                                plt.text(x[0],(y[1]*(0.1)),a[0],rotation=0,color='b',fontsize=20,horizontalalignment='right',fontweight='bold')
                                        else:
                                                plt.text(x[0],(y[1]*(0.1)),a[0],rotation=0,color='r',fontsize=20,horizontalalignment='right',fontweight='bold')
                                elif len(a[0])<3:
                                        plt.plot(x,y,c,linewidth=0.3,alpha=0.2,linestyle='-.')
                                        if a[0][-1].upper()=='S':
                                                plt.text(x[0],(y[1]*(0.4)) * (-1),a[0],rotation=-0,color='b',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',alpha=0.5)
                                        else:
                                                plt.text(x[0],(y[1]*(0.4)) * (1),a[0],rotation=-0,color='r',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',alpha=0.5)
                                elif len(a[0])==3:
                                        plt.plot(x,y,c,linewidth=0.3,alpha=0.2,linestyle='-.')
                                        if a[0][-1].upper()=='S':
                                                plt.text(x[0],(y[1]*(0.5)) * (-1),a[0],rotation=-0,color='b',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',alpha=0.5,verticalalignment='bottom')
                                        else:
                                                plt.text(x[0],(y[1]*(0.5)) * (1),a[0],rotation=-0,color='r',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',alpha=0.5,verticalalignment='top')
                                elif (len([m.start() for m in re.finditer('I',a[0].upper())])>=1):
                                        plt.plot(x,y,'gray',linewidth=0.3,alpha=0.9,linestyle='-.')
                                        if a[0][-1].upper()=='S':
                                                plt.text(x[0],(y[1]*(1)) * (-1),a[0],rotation=90,color='gray',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',verticalalignment='bottom',alpha=0.5)
                                        else:
                                                plt.text(x[0],(y[1]*(1)) * (1),a[0],rotation=90,color='gray',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',verticalalignment='top',alpha=0.5)
                                elif (len([m.start() for m in re.finditer('K',a[0].upper())])>=1):
                                        plt.plot(x,y,'gray',linewidth=0.3,alpha=0.9,linestyle='-.')
                                        if a[0][-1].upper()=='S':
                                                plt.text(x[0],(y[1]*(0.7)) * (-1),a[0],rotation=90,color='gray',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',verticalalignment='bottom')
                                        else:
                                                plt.text(x[0],(y[1]*(0.6)) * (1),a[0],rotation=90,color='gray',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',verticalalignment='top')

def map(color_by=None,stations_subset=None,stations=None,figsize=None,key=0,legend_title='Legend',legx='"auto"',legy='"auto"',hexclrs=None):
        sta_set = None
        sta_subset = None
        if stations_subset is None and stations is None:
                sta_set,sta_subset = ObsQA.io.getstalist()
        if stations is not None:
                sta_set = stations.copy()
        if stations_subset is not None:
                sta_subset = stations_subset.copy()
        if sta_set is not None:
                labnds = [sta_set['Latitude'].min(),sta_set['Latitude'].max()]
                lobnds = [sta_set['Longitude'].min(),sta_set['Longitude'].max()]
        if sta_subset is not None:
                labnds = [sta_subset['Latitude'].min(),sta_subset['Latitude'].max()]
                lobnds = [sta_subset['Longitude'].min(),sta_subset['Longitude'].max()]
        sta_set_orig = sta_set.copy()
        if color_by is not None:
                items = sta_set[color_by].unique()
                sta_set = [sta_set.iloc[np.in1d(sta_set[color_by],item)] for item in items]
        else:
                sta_set = [sta_set]
        tilekey = dict()
        # good map tiles: http://leaflet-extras.github.io/leaflet-providers/preview/
        tilekey[0] = ('https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}','Tiles &copy; Esri &mdash; Sources: GEBCO, NOAA, CHS, OSU, UNH, CSUMB, National Geographic')
        tilekey[1] = ('https://server.arcgisonline.com/ArcGIS/rest/services/World_Terrain_Base/MapServer/tile/{z}/{y}/{x}','Tiles &copy; Esri &mdash; Source: USGS, Esri, TANA, DeLorme, and NPS')
        tilekey[2] = ('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}','Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community')
        if figsize is None:
                if sta_subset is not None:
                        figsize = [600,500]
                else:
                        figsize = [2000,800]
        fig = Figure(width=figsize[0], height=figsize[1])
        m = folium.Map(tiles=tilekey[key][0],attr=tilekey[key][1],zoom_start=1,min_zoom=1,minlon = lobnds[0],maxlon=lobnds[1],control=False,max_bounds=True,zoom_control=False,scrollWheelZoom=True) #,width=figsize[0],height=figsize[1]
        m.add_child(folium.LatLngPopup())
        if hexclrs is None:
                if len(sta_set)==1:
                        hexclrs = ['#e41a1c'] #red
                else:
                        colors = distinctipy.get_colors(len(sta_set),pastel_factor=0.0,exclude_colors=[(1.0, 0.0, 1.0),(1,1,1),(0,0,0),(0,1,0),(0,1,1)],colorblind_type='Deuteranomaly',n_attempts=1000)
                        hexclrs = [matplotlib.colors.to_hex(colors[i]) for i in range(len(colors))][0:len(colors)]
        hexclrs = hexclrs[0:len(sta_set)]
        if color_by is None:
                color_by = 'Station'
        if sta_set is not None:
                legend_dict = dict()
                for clri,cur_sta_set in enumerate(sta_set):
                        color = hexclrs[clri]
                        label = str(cur_sta_set[color_by].tolist()[0])
                        legend_dict[label]  = color
                        locationlist = (np.array((cur_sta_set['Latitude'],list(cur_sta_set['Longitude']))).T.tolist())
                        g1 = folium.FeatureGroup(name=('<span style=>{txt}</span>').format( txt='All Stations'))
                        for point in range(0, len(locationlist)):
                                marker = folium.Marker(locationlist[point], popup=list(cur_sta_set['Network'] + '.' + list(cur_sta_set['Station']) + ' (' + list(cur_sta_set['Experiment']) + ')')[point],icon=folium.Icon(color='white',icon_color=color,prefix='fa',icon='caret-up',shadow_size=0,icon_size=(0,0),icon_anchor=(0,0)),legend_name='hello')
                                L = folium.FeatureGroup().add_child(marker)
                                marker.add_to(g1)
                        m.add_child(g1)
                if len(sta_set)==1:
                        labels = [color_by]
                else:
                        labels = list(legend_dict.keys())
                colors = list(legend_dict.values())
                macro = ObsQA.getLegendTemplate(colors,labels,title_str=legend_title,legx=legx,legy=legy)
                m.add_child(macro)
        if sta_subset is not None:
                color = 'red'
                locationlist_subset = (np.array((sta_subset['Latitude'],list(sta_subset['Longitude']))).T.tolist())
                g2 = folium.FeatureGroup(name=('<span style=>{txt}</span>').format( txt='Subset'),keep_in_front=True)
                for point in range(0, len(locationlist_subset)):
                        marker = folium.Marker(locationlist_subset[point], popup=list(sta_subset['Network'] + '.' + list(sta_subset['Station']) + ' (' + list(sta_subset['Experiment']) + ')')[point],icon=folium.Icon(color='white',icon_color=color,prefix='fa',icon='caret-up',shadow_size=0,icon_size=(0,0),icon_anchor=(0,0)),keep_in_front=True)
                        marker.add_to(g2)
                m.add_child(g2)
        if sta_set is not None:
                fit_bounds=[sta_set_orig[['Latitude', 'Longitude']].min().values.tolist(), sta_set_orig[['Latitude', 'Longitude']].max().values.tolist()]
                m.fit_bounds(fit_bounds)
        if sta_subset is not None:
                fit_bounds_subset = [sta_subset[['Latitude','Longitude']].min().tolist(), sta_subset[['Latitude','Longitude']].max().tolist()]
                m.fit_bounds(fit_bounds_subset)
        # if sta_set is not None:
        #         m.add_child(macro)
        fig.add_child(m)
        return m

def mapby(df,color_by=None,height=900,width=1000,map_set='NOAA/NGDC/ETOPO1',elevation=['black','white', 'black','white'],layercontrol=False,center=None,file=None):
        if color_by is not None:
                items = df[color_by].unique()
                sta_set = [df.iloc[np.in1d(df[color_by],item)] for item in items]
        else:
                sta_set = [df]
        MapObj = geemap.Map(height=height,width=width)
        landcover = ee.Image(map_set).select('bedrock')
        elevationVis = {'min': -7000.0,'max': 3000.0,'palette': elevation}
        MapObj.addLayer(landcover, elevationVis, 'NOAA')
        if layercontrol:
                # MapObj.addLayerControl()
                pass
        if len(sta_set)==1:
                hexclrs = ['#e41a1c'] #red
        else:
                colors = distinctipy.get_colors(len(sta_set),pastel_factor=0.0,exclude_colors=[(1.0, 0.0, 1.0),(1,1,1),(0,0,0),(0,1,0),(0,1,1)],colorblind_type='Deuteranomaly',n_attempts=1000)
                hexclrs = [matplotlib.colors.to_hex(colors[i]) for i in range(len(colors))]
        legend_dict = dict()
        for clri,cur_sta_set in enumerate(sta_set):
                color = hexclrs[clri]
                if color_by is not None:
                        label = cur_sta_set[color_by].tolist()[0]
                        legend_dict[label]  = color
                        MapObj.add_circle_markers_from_xy(cur_sta_set, x="Longitude", y="Latitude", radius=1, color=color, fill_color="black",label=label)
                else:
                        MapObj.add_circle_markers_from_xy(cur_sta_set, x="Longitude", y="Latitude", radius=1, color=color, fill_color="black")
                        # MapObj.add_legend(legend_title="Test", legend_dict=legend_dict,position="bottomright")
        if color_by is not None:
                MapObj.add_legend(title=color_by, legend_dict=legend_dict,position="bottomright")
        if center is not None:
                MapObj.set_center(lat=center[0],lon=center[1],zoom=center[2])
        else:
                MapObj.zoom_to_bounds((df.Longitude.min(),df.Latitude.min(),df.Longitude.max(),df.Latitude.max()))
        MapObj.set_control_visibility(layerControl=False,fullscreenControl=False,latLngPopup=False)
        if file is not None:
                MapObj
                MapObj.to_image(file)
        return MapObj

def ML_fig_event_corrected(PreProcEvent,CorrectedEvent,evstream=None, TF_list=None, fmin=1./150., fmax=2.,prefix = '',yes_filter=True,ylon=True,scale='linear',yhard=None,spl_list = [611,612,613,614,615,616]):
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
        static_order = ['Z1','Z2-1','ZP-21','ZH','ZP-H','ZP']
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
        # plt.figure(figsize=(8, 8))
        eventtime = CorrectedEvent.eventid.squeeze().to_list()[0]
        evtstamp = str(UTCDateTime(eventtime).year) + '.' + str(UTCDateTime(eventtime).julday) + '.' + str(UTCDateTime(eventtime).hour) + '.' + str(UTCDateTime(eventtime).minute)
        i = 0
        for key in static_order:
                plt.subplot(int(str(spl_list[i])[0]),int(str(spl_list[i])[1]),int(str(spl_list[i])[2:]))
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

def Py_fig_event_corrected(evstream, TF_list, fmin=1./150., fmax=2.,prefix='',spl_list = [611,612,613,614,615,616]):
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
        static_order = ['Z1','Z2-1','ZP-21','ZH','ZP-H','ZP']
        i = 0
        for key in static_order:
                plt.subplot(int(str(spl_list[i])[0]),int(str(spl_list[i])[1]),int(str(spl_list[i])[2:]))
                if TF_list[key]:
                        plt.plot(taxis, trZ.data, 'lightgray', lw=0.5)
                        tr = Trace(data=evstream.correct[key],header=trZ.stats).filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
                        plt.plot(taxis, tr.data, 'k', lw=0.5)
                plt.title(prefix + evstream.key + ' ' + evstream.tstamp + ': ' + key, fontdict={'fontsize': 8})
                plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,scilimits=(-3, 3))
                plt.xlim((0., trZ.stats.npts/sr))
                i +=1
        plt.xlabel('Time since earthquake (sec)')
        plt.tight_layout()

def TraceResiduals(event_time,network,sta,py_folders,ml_folders,respct=False,scale='linear',spl_list = [611,612,613,614,615,616]):
        '''Just writes residuals to the obspy trace objects to feed into ML_fig_event_corrected'''
        eventtime_jdaystr = str(UTCDateTime(event_time).year) + '.' + str(UTCDateTime(event_time).julday).zfill(3) + '.' + str(UTCDateTime(event_time).hour).zfill(2) + '.' + str(UTCDateTime(event_time).minute).zfill(2)
        PreProcEvent,CorrectedEvent,ML_TFs = ObsQA.io.GetML_EventData_and_TransferFunctions(event_time,network,sta,ml_folders[0],ml_folders[1],ml_folders[2])
        for i in range(len(PreProcEvent)):
                tmp = PreProcEvent.iloc[i]['Trace']
                tmp.data = tmp.data*0
                PreProcEvent.iat[i,-1] = tmp
        ML_CorrectedEvent = CorrectedEvent
        path = py_folders[1] + '/' + network + '.' + sta + '/CORRECTED/' + network + '.' + sta + '.' + eventtime_jdaystr + '.day'
        py_files = g.glob(path + '*.pkl')
        f = py_files[0]
        evstream = pkl.load(open(f,'rb'))
        keys = CorrectedEvent['label'].to_list()
        delta = pd.DataFrame.from_dict({'label':keys,'residual':[CorrectedEvent.iloc[0].Trace.data*0 for i in range(len(keys))]})
        for k in keys:
                PyTr = evstream.correct[k]
                MLTr = CorrectedEvent.iloc[(CorrectedEvent['label']==k).to_list()].timeseries.to_list()[0]
                if respct:
                        res = ((PyTr - MLTr)/MLTr)
                else:
                        res = PyTr - MLTr
                i = np.where(delta['label']==k)[0][0]
                delta.at[i,'residual'] = res
                tmp = CorrectedEvent.iloc[i,-1]
                tmp.data = delta.at[i,'residual']
                CorrectedEvent.iloc[i,-1] = tmp
        prefix = '<Py-ML Residual Fraction>'
        TF_list = {i : j for i, j in zip(CorrectedEvent.label.to_list(), np.ones(len(CorrectedEvent.label.to_list()),dtype=bool))}
        ObsQA.plots.ML_fig_event_corrected(PreProcEvent,CorrectedEvent,evstream=None, TF_list=TF_list, fmin=1./150., fmax=2.,prefix=prefix,yes_filter=False,ylon=False,scale=scale,spl_list=spl_list)

def Comp_PyML_EventData(event_time,network,stalist,py_folders,ml_folders):
        spl_list = [631,634,637,6310,6313,6316]
        for stai in stalist:
                # Matlab ATaCR output plot
                prefix = '<Matlab-ATaCR-Output> '
                PreProcEvent,CorrectedEvent,ML_TFs = ObsQA.io.GetML_EventData_and_TransferFunctions(event_time,network,stai,ml_folders[0],ml_folders[1],ml_folders[2])
                if stai=='M02A':
                        ylon = False
                        yhard = [-0.2*np.mean(np.array(PreProcEvent['data'].to_list())),0.2*np.mean(PreProcEvent['data'].to_list())]
                else:
                        ylon = True
                        yhard = None
                TF_list = {i : j for i, j in zip(CorrectedEvent.label.to_list(), np.ones(len(CorrectedEvent.label.to_list()),dtype=bool))}
                plt.figure(figsize=(18, 12))
                ObsQA.plots.ML_fig_event_corrected(PreProcEvent,CorrectedEvent,evstream=None, TF_list=TF_list, fmin=1./150., fmax=2.,prefix=prefix,ylon = ylon,yhard = yhard,spl_list=list(np.array(spl_list)+0))
                # Python ATaCR output plot
                eventtime_jdaystr = str(UTCDateTime(event_time).year) + '.' + str(UTCDateTime(event_time).julday).zfill(3) + '.' + str(UTCDateTime(event_time).hour).zfill(2) + '.' + str(UTCDateTime(event_time).minute).zfill(2)
                path = py_folders[1] + '/' + network + '.' + stai + '/CORRECTED/' + network + '.' + stai + '.' + eventtime_jdaystr + '.day'
                py_files = g.glob(path + '*.pkl')
                f = py_files[0]
                display('Matlab File: ' + PreProcEvent.iloc[0].Folder + PreProcEvent.iloc[0].File)
                display('Python File : ' + f)
                evstream = pkl.load(open(f,'rb'))
                TF_list = {i : j for i, j in zip(list(evstream.correct.keys()), np.ones(len(list(evstream.correct.keys())),dtype=bool))}
                prefix = '<Python-ATaCR-Output> '
                ObsQA.plots.Py_fig_event_corrected(evstream, TF_list, prefix=prefix,spl_list=list(np.array(spl_list)+1))
                # Residual plot
                ObsQA.plots.TraceResiduals(event_time,network,stai,py_folders,ml_folders,respct=True,scale='symlog',spl_list=list(np.array(spl_list)+2))
                plt.show()
                print('---------------------------------------------------------------------------------------------------------')

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
                        ph = ObsQA.qa.Phase(M.__dict__[a])
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
                                if i==(len(ph_map)-1):
                                        plt.xlabel('Frequency (Hz)')
                                if i<(len(ph_map)-1):
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
                        ad = ObsQA.qa.Admittance(M.__dict__[a],M.__dict__[b])
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
                        coh = ObsQA.qa.Coherence(M.__dict__[a],M.__dict__[b],M.__dict__[c])
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