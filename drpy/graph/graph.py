import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cmx
import matplotlib.colors as colors 
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.patheffects as PathEffects
import cartopy.io.shapereader as shpreader
import pandas as pd 
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.colors as colors
import matplotlib.patheffects as PathEffects
from . import colormaps_drpy as cmaps

def make_colorbar(ax,vmin,vmax,cmap):
    cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap,
                              norm=cNorm,
                              orientation='horizontal',extend='both')
    return cb1

class case_study:

  def __init__(self,filename=None,center_lat=None,center_lon=None,path_to_models='../models/'):

    import drpy 
    dpr = drpy.core.GPMDPR(filename=filename)
    dpr.read()
    dpr.parse_dtime()
    #if no center point is given, use middle of orbit. 
    if (center_lat is None) or (center_lon is None):
        #determine map center
        s = 0 
        e = dpr.ds.Longitude.shape[0]
        middle = int((e-s)/2)
        center_lon = dpr.ds.Longitude.values[middle,24]
        center_lat = dpr.ds.Latitude.values[middle,24]
    corners = [center_lon - 5,center_lon + 5,center_lat-5,center_lat+5]
    dpr.corners = corners
    dpr.setboxcoords()
    #drop dead weight (i.e. blank data)
    dpr.ds = dpr.ds.dropna(dim='nscan',how='all')
    self.dpr = dpr
    self.path_to_models=path_to_models

  def plotter_along(self,start_index=25,end_index=-25,scan=24,params_new=None):

    if params_new is None:
      params = {'z_vmin':10,'z_vmax':40,'y_max':10,'dfr_vmin':-2,'dfr_vmax':10,
                'dm_vmin':0,'dm_vmax':2,'nw_vmin':1,'nw_vmax':6,'r_vmin':-1,'r_vmax':2,
                'xsections':[0,2,4],'temperature':False,'t_levels':np.linspace(-20,20,10)}
    else:
      params = {'z_vmin':10,'z_vmax':40,'y_max':10,'dfr_vmin':-2,'dfr_vmax':10,
                'dm_vmin':0,'dm_vmax':2,'nw_vmin':1,'nw_vmax':6,'r_vmin':-1,'r_vmax':2,
                'xsections':[0,2,4],'temperature':False,'t_levels':np.linspace(-20,20,10)}
      keys_old = list(params.keys())
      keys = list(params_new.keys())
      for key in keys: 
        for i in keys_old:
                if key in i:
                  params[key] = params_new[key]


    import proplot as plot
    #plot parameters that I personally like, feel free to make these your own.
    matplotlib.rcParams['axes.facecolor'] = [0.9,0.9,0.9]
    matplotlib.rcParams['axes.labelsize'] = 14
    matplotlib.rcParams['axes.titlesize'] = 14
    matplotlib.rcParams['xtick.labelsize'] = 12
    matplotlib.rcParams['ytick.labelsize'] = 12
    matplotlib.rcParams['legend.fontsize'] = 12
    matplotlib.rcParams['legend.facecolor'] = 'w'
    matplotlib.rcParams['savefig.transparent'] = False

    #determine map center
    s = 0 
    e = self.dpr.ds.Longitude.shape[0]
    middle = int((e-s)/2)
    lon0 = self.dpr.ds.Longitude.values[middle,24]
    lat0 = self.dpr.ds.Latitude.values[middle,24]
    #set specific aspect for cartopy plot
    x = 7.5
    y = 0.6666666*x
    #set corners 
    corners = [lon0-x,lon0+x,lat0-y,lat0+y]

    #draw axes with proplot
    array = [  # the "picture" (1 == subplot A, 2 == subplot B, etc.)
    [1,1,1,0,0,0,0],
    [1,1,1,4,4,4,4],
    [2,2,2,4,4,4,4],
    [2,2,2,4,4,4,4],
    [3,3,3,4,4,4,4],
    [3,3,3,0,0,0,0]]

    fig, axs = plot.subplots(array, width=10,height=5,span=False,proj=['cart', 'cart', 'cart','cyl',],tight=False)

    fig.set_facecolor('w')

    #do things to map subplot 
    ax = axs[3]
    ax.format(borderscolor='k',gridminor=True,
              borders=True,lonlabels='t', latlabels='l',
              lonlim=(corners[0], corners[1]), latlim=(corners[2], corners[3]))
    
    #manually adjust map location to make it fit nicely 
    box = ax.get_position()
    box.y0 = box.y0 - 0.05
    box.y1 = box.y1 - 0.05
    box.x0 = box.x0 + 0.01
    box.x1 = box.x1 + 0.01
    ax.set_position(box)
    # add land and ocean colors
    ax.add_feature(cartopy.feature.LAND,facecolor=[0.9,0.9,0.9])
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.STATES)
    plt.setp(ax.spines.values(), color='orangered',lw=2)

    #!!!ADD ADDITIONAL MAP STUFF HERE!!!# 
    
    #plot swath 
    ax.plot(self.dpr.ds.Longitude[:,0]+0.0485,self.dpr.ds.Latitude[:,0],'--k')
    ax.plot(self.dpr.ds.Longitude[:,-1]-0.0485,self.dpr.ds.Latitude[:,-1],'--k')

    #plot data on map 
    pm = ax.scatter(self.dpr.ds.Longitude[:,:],self.dpr.ds.Latitude[:,:],c=self.dpr.ds.zFactorFinalNearSurface[:,:,0],vmin=params['z_vmin'],vmax=params['z_vmax'],s=1,cmap='Spectral_r',linewidths=0.1,zorder=10)
    s = start_index
    e = end_index
    w = scan

    #plot the along_track choice
    ax.plot(self.dpr.ds.Longitude[s:e,w],self.dpr.ds.Latitude[s:e,w],'-k',markerfacecolor='w',ms=10,zorder=12)
    ax.plot(self.dpr.ds.Longitude[s,w],self.dpr.ds.Latitude[s,w],'k',markerfacecolor='w',ms=10,label='Start',marker='$L$',markeredgewidth=1,zorder=12)
    ax.plot(self.dpr.ds.Longitude[e,w],self.dpr.ds.Latitude[e,w],'k',markerfacecolor='w',ms=10,label='End',marker='$R$',markeredgewidth=1,zorder=12)

    #add coastlines 
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'),color='k',zorder=11)

    #add zoomed out map for context using proplot
    inset_axis = ax.inset([0.575,-0.4,0.5,0.5],proj='cyl',zoom=False,zorder=11)
    inset_axis.format(land=True,landcolor=[0.9,0.9,0.9],borderscolor='k',gridminor=True)
    inset_axis.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
    inset_axis.plot(self.dpr.ds.Longitude[:,0]+0.0485,self.dpr.ds.Latitude[:,0],'--k',lw=0.5,)
    inset_axis.plot(self.dpr.ds.Longitude[:,-1]-0.0485,self.dpr.ds.Latitude[:,-1],'--k',lw=0.5,)
    inset_axis.plot([corners[0],corners[0],corners[1],corners[1],corners[0]],[corners[2],corners[3],corners[3],corners[2],corners[2]],'-',color='orangered')
    timestr = pd.to_datetime(self.dpr.ds.time[middle,24].values).strftime(format='%Y-%m-%d %H:%M')
    text = inset_axis.text(-0.05,-0.2,'Scan Time: ' + timestr,transform=ax.transAxes,fontsize=10)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
    text = inset_axis.text(0.025,-0.275,'Created with DRpy',transform=ax.transAxes,fontsize=10)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])

    #draw colorbars in better spot (above map)
    ax_cbar1 = fig.add_axes([0.6, 0.93, 0.33, 0.015])
    ax_cbar2 = fig.add_axes([0.6, 0.85, 0.33, 0.015])
    cbar1=True
    cbar2=True
    zfilled=True
    dfrfilled=True
    dmfilled=True
    nwfilled=True
    rfilled=True
    for x in np.arange(0,3):
        xsection_i = params['xsections'][x]
        if (xsection_i <=3) and (cbar1 or cbar2) and (zfilled):
            if cbar1:
                ax_cbar = ax_cbar1
                cbar1=False
            else:
                ax_cbar = ax_cbar2
                cbar2=False
            #add Z colorbar 
            text = ax_cbar.text(-0.25,0,'$Z_{e}$, [$dBZ$]',fontsize=10,transform=ax_cbar.transAxes)
            text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
            cb1 = make_colorbar(ax_cbar,params['z_vmin'],params['z_vmax'],plt.cm.Spectral_r)
            zfilled=False
        elif ((xsection_i==4) or (xsection_i==5)) and (cbar1 or cbar2) and (dfrfilled):
            if cbar1:
                ax_cbar = ax_cbar1
                cbar1=False
            else:
                ax_cbar = ax_cbar2
                cbar2=False
            text = ax_cbar.text(-0.33,0,'$DFR_{Ku-Ka}$, [$dB$]',fontsize=10,transform=ax_cbar.transAxes)
            text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
            cb2 = make_colorbar(ax_cbar,params['dfr_vmin'],params['dfr_vmax'],cmaps.turbo)
            dfrfilled=False
        elif (xsection_i==6) or (xsection_i==9) or (xsection_i==10) and (cbar1 or cbar2) and (dmfilled):
            if cbar1:
                ax_cbar = ax_cbar1
                cbar1=False
            else:
                ax_cbar = ax_cbar2
                cbar2=False
            text = ax_cbar.text(-0.33,0,'$D_m$, [$mm$]',fontsize=10,transform=ax_cbar.transAxes)
            text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
            cb2 = make_colorbar(ax_cbar,params['dm_vmin'],params['dm_vmax'],cmaps.plasma)
            dmfilled=False
        elif (xsection_i==7) and (cbar1 or cbar2) and (nwfilled):
            if cbar1:
                ax_cbar = ax_cbar1
                cbar1=False
            else:
                ax_cbar = ax_cbar2
                cbar2=False
            text = ax_cbar.text(-0.33,0,'$N_w$, [$\log{(mm^{-1} \ m^{-3})}$]',fontsize=10,transform=ax_cbar.transAxes)
            text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
            cb2 = make_colorbar(ax_cbar,params['nw_vmin'],params['nw_vmax'],cmaps.plasma)
            nwfilled=False
        elif (xsection_i==8) or (xsection_i==11) and (cbar1 or cbar2) and (rfilled):
            if cbar1:
                ax_cbar = ax_cbar1
                cbar1=False
            else:
                ax_cbar = ax_cbar2
                cbar2=False
            text = ax_cbar.text(-0.33,0,'$R$, [$\log{(mm \ hr^{-1})}$]',fontsize=10,transform=ax_cbar.transAxes)
            text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
            cb2 = make_colorbar(ax_cbar,params['r_vmin'],params['r_vmax'],cmaps.plasma)
            rfilled=False


    #get distances to plot cross-section 
    self.dpr.get_physcial_distance(reference_point=[self.dpr.ds.Longitude.values[s,w],self.dpr.ds.Latitude.values[s,w]])

    if np.max(params['xsections']) >= 9:
        print('running Chase retrieval. Please wait...')
        self.dpr.run_Chase2021(models_path=self.path_to_models)
        print('done')

    for i in np.arange(0,3):
      ax = axs[i]
      cross_section_method(self,params['xsections'][i],s,e,w,ax,params)

    axs[1].set_ylabel('Altitude, [km]',labelpad=8)
    ax.set_xlabel('Along Track Distance, [km]',labelpad=8)

  def plotter_cross(self,along_track_index=25,params_new=None):

    if params_new is None:
      params = {'z_vmin':10,'z_vmax':40,'y_max':10,'dfr_vmin':-2,'dfr_vmax':10,
                'dm_vmin':0,'dm_vmax':2,'nw_vmin':1,'nw_vmax':6,'r_vmin':-1,'r_vmax':2,
                'xsections':[0,2,4],'temperature':False,'t_levels':np.linspace(-20,20,10)}
    else:
      params = {'z_vmin':10,'z_vmax':40,'y_max':10,'dfr_vmin':-2,'dfr_vmax':10,
                'dm_vmin':0,'dm_vmax':2,'nw_vmin':1,'nw_vmax':6,'r_vmin':-1,'r_vmax':2,
                'xsections':[0,2,4],'temperature':False,'t_levels':np.linspace(-20,20,10)}
      keys_old = list(params.keys())
      keys = list(params_new.keys())
      for key in keys: 
        for i in keys_old:
                if key in i:
                  params[key] = params_new[key]


    import proplot as plot
    #plot parameters that I personally like, feel free to make these your own.
    matplotlib.rcParams['axes.facecolor'] = [0.9,0.9,0.9]
    matplotlib.rcParams['axes.labelsize'] = 14
    matplotlib.rcParams['axes.titlesize'] = 14
    matplotlib.rcParams['xtick.labelsize'] = 12
    matplotlib.rcParams['ytick.labelsize'] = 12
    matplotlib.rcParams['legend.fontsize'] = 12
    matplotlib.rcParams['legend.facecolor'] = 'w'
    matplotlib.rcParams['savefig.transparent'] = False

    #determine map center
    s = 0 
    e = self.dpr.ds.Longitude.shape[0]
    middle = int((e-s)/2)
    lon0 = self.dpr.ds.Longitude.values[middle,24]
    lat0 = self.dpr.ds.Latitude.values[middle,24]
    #set specific aspect for cartopy plot
    x = 7.5
    y = 0.6666666*x
    #set corners 
    corners = [lon0-x,lon0+x,lat0-y,lat0+y]

    #draw axes with proplot
    array = [  # the "picture" (1 == subplot A, 2 == subplot B, etc.)
    [1,1,1,0,0,0,0],
    [1,1,1,4,4,4,4],
    [2,2,2,4,4,4,4],
    [2,2,2,4,4,4,4],
    [3,3,3,4,4,4,4],
    [3,3,3,0,0,0,0]]

    fig, axs = plot.subplots(array, width=10,height=5,span=False,proj=['cart', 'cart', 'cart','cyl',],tight=False)

    fig.set_facecolor('w')

    #do things to map subplot 
    ax = axs[3]
    ax.format(borderscolor='k',gridminor=True,
              borders=True,lonlabels='t', latlabels='l',
              lonlim=(corners[0], corners[1]), latlim=(corners[2], corners[3]))
    
    #manually adjust map location to make it fit nicely 
    box = ax.get_position()
    box.y0 = box.y0 - 0.05
    box.y1 = box.y1 - 0.05
    box.x0 = box.x0 + 0.01
    box.x1 = box.x1 + 0.01
    ax.set_position(box)
    # add land and ocean colors
    ax.add_feature(cartopy.feature.LAND,facecolor=[0.9,0.9,0.9])
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.STATES)
    plt.setp(ax.spines.values(), color='orangered',lw=2)

    #!!!ADD ADDITIONAL MAP STUFF HERE!!!# 
    
    #plot swath 
    ax.plot(self.dpr.ds.Longitude[:,0]+0.0485,self.dpr.ds.Latitude[:,0],'--k')
    ax.plot(self.dpr.ds.Longitude[:,-1]-0.0485,self.dpr.ds.Latitude[:,-1],'--k')

    #plot data on map 
    pm = ax.scatter(self.dpr.ds.Longitude[:,:],self.dpr.ds.Latitude[:,:],c=self.dpr.ds.zFactorFinalNearSurface[:,:,0],vmin=params['z_vmin'],vmax=params['z_vmax'],s=1,cmap='Spectral_r',linewidths=0.1,zorder=10)
    a = along_track_index

    #plot the along_track choice
    ax.plot(self.dpr.ds.Longitude[a,:],self.dpr.ds.Latitude[a,:],'-k',markerfacecolor='w',ms=10,zorder=12)
    ax.plot(self.dpr.ds.Longitude[a,0],self.dpr.ds.Latitude[a,0],'k',markerfacecolor='w',ms=10,label='Start',marker='$L$',markeredgewidth=1,zorder=12)
    ax.plot(self.dpr.ds.Longitude[a,-1],self.dpr.ds.Latitude[a,-1],'k',markerfacecolor='w',ms=10,label='End',marker='$R$',markeredgewidth=1,zorder=12)

    #add coastlines 
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'),color='k',zorder=11)

    #add zoomed out map for context using proplot
    inset_axis = ax.inset([0.575,-0.4,0.5,0.5],proj='cyl',zoom=False,zorder=11)
    inset_axis.format(land=True,landcolor=[0.9,0.9,0.9],borderscolor='k',gridminor=True)
    inset_axis.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
    inset_axis.plot(self.dpr.ds.Longitude[:,0]+0.0485,self.dpr.ds.Latitude[:,0],'--k',lw=0.5,)
    inset_axis.plot(self.dpr.ds.Longitude[:,-1]-0.0485,self.dpr.ds.Latitude[:,-1],'--k',lw=0.5,)
    inset_axis.plot([corners[0],corners[0],corners[1],corners[1],corners[0]],[corners[2],corners[3],corners[3],corners[2],corners[2]],'-',color='orangered')
    timestr = pd.to_datetime(self.dpr.ds.time[middle,24].values).strftime(format='%Y-%m-%d %H:%M')
    text = inset_axis.text(-0.05,-0.2,'Scan Time: ' + timestr,transform=ax.transAxes,fontsize=10)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
    text = inset_axis.text(0.025,-0.275,'Created with DRpy',transform=ax.transAxes,fontsize=10)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])

    #draw colorbars in better spot (above map)
    ax_cbar1 = fig.add_axes([0.6, 0.93, 0.33, 0.015])
    ax_cbar2 = fig.add_axes([0.6, 0.85, 0.33, 0.015])
    cbar1=True
    cbar2=True
    zfilled=True
    dfrfilled=True
    dmfilled=True
    nwfilled=True
    rfilled=True
    for x in np.arange(0,3):
        xsection_i = params['xsections'][x]
        if (xsection_i <=3) and (cbar1 or cbar2) and (zfilled):
            if cbar1:
                ax_cbar = ax_cbar1
                cbar1=False
            else:
                ax_cbar = ax_cbar2
                cbar2=False
            #add Z colorbar 
            text = ax_cbar.text(-0.25,0,'$Z_{e}$, [$dBZ$]',fontsize=10,transform=ax_cbar.transAxes)
            text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
            cb1 = make_colorbar(ax_cbar,params['z_vmin'],params['z_vmax'],plt.cm.Spectral_r)
            zfilled=False
        elif ((xsection_i==4) or (xsection_i==5)) and (cbar1 or cbar2) and (dfrfilled):
            if cbar1:
                ax_cbar = ax_cbar1
                cbar1=False
            else:
                ax_cbar = ax_cbar2
                cbar2=False
            text = ax_cbar.text(-0.33,0,'$DFR_{Ku-Ka}$, [$dB$]',fontsize=10,transform=ax_cbar.transAxes)
            text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
            cb2 = make_colorbar(ax_cbar,params['dfr_vmin'],params['dfr_vmax'],cmaps.turbo)
            dfrfilled=False
        elif (xsection_i==6) or (xsection_i==9) or (xsection_i==10) and (cbar1 or cbar2) and (dmfilled):
            if cbar1:
                ax_cbar = ax_cbar1
                cbar1=False
            else:
                ax_cbar = ax_cbar2
                cbar2=False
            text = ax_cbar.text(-0.33,0,'$D_m$, [$mm$]',fontsize=10,transform=ax_cbar.transAxes)
            text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
            cb2 = make_colorbar(ax_cbar,params['dm_vmin'],params['dm_vmax'],cmaps.plasma)
            dmfilled=False
        elif (xsection_i==7) and (cbar1 or cbar2) and (nwfilled):
            if cbar1:
                ax_cbar = ax_cbar1
                cbar1=False
            else:
                ax_cbar = ax_cbar2
                cbar2=False
            text = ax_cbar.text(-0.33,0,'$N_w$, [$\log{(mm^{-1} \ m^{-3})}$]',fontsize=10,transform=ax_cbar.transAxes)
            text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
            cb2 = make_colorbar(ax_cbar,params['nw_vmin'],params['nw_vmax'],cmaps.plasma)
            nwfilled=False
        elif (xsection_i==8) or (xsection_i==11) and (cbar1 or cbar2) and (rfilled):
            if cbar1:
                ax_cbar = ax_cbar1
                cbar1=False
            else:
                ax_cbar = ax_cbar2
                cbar2=False
            text = ax_cbar.text(-0.33,0,'$R$, [$\log{(mm \ hr^{-1})}$]',fontsize=10,transform=ax_cbar.transAxes)
            text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
            cb2 = make_colorbar(ax_cbar,params['r_vmin'],params['r_vmax'],cmaps.plasma)
            rfilled=False


    #get distances to plot cross-section 
    self.dpr.get_physcial_distance(reference_point=[self.dpr.ds.Longitude.values[a,0],self.dpr.ds.Latitude.values[a,0]])

    if np.max(params['xsections']) >= 9:
        print('running Chase retrieval. Please wait...')
        self.dpr.run_Chase2021(models_path=self.path_to_models)
        print('done')

    for i in np.arange(0,3):
      ax = axs[i]
      cross_section_method2(self,params['xsections'][i],along_track_index,ax,params)

    axs[1].set_ylabel('Altitude, [km]',labelpad=8)
    ax.set_xlabel('Cross Track Distance, [km]',labelpad=8)
#xsection choice key
# 0   raw Ku band 
# 1   corrected Ku band 
# 2   raw Ka band 
# 3   corrected Ka band 
# 4   raw DFR        (rawKu - rawKa)
# 5   corrected DFR  (Ku - Ka)
# 6   retrieved Dm_l (2A.DPR method)
# 7   retrieved Nw   (2A.DPR method)
# 8   retreived R    (2A.DPR method)
# 9   retrieved Dm_l (Chase method)
# 10  retrieved Dm_s (Chase method)
# 11  retrieved R    (Chase method)


def cross_section_method(self,choice,s,e,w,ax,params): 
    from drpy.graph import cmaps
    if choice ==0:
        x = self.dpr.ds.zFactorMeasured[:,:,:,0].where(self.dpr.ds.zFactorMeasured[:,:,:,0] >= 10)
        cmap = 'Spectral_r'
        levels=np.linspace(params['z_vmin'], params['z_vmax'], 50)
        vmin = params['z_vmin']
        vmax = params['z_vmax']
        label = 'KuPR'
    elif choice==1:
        x = self.dpr.ds.zFactorFinal[:,:,:,0].where(self.dpr.ds.zFactorFinal[:,:,:,0] >= 10)
        cmap = 'Spectral_r'
        levels=np.linspace(params['z_vmin'], params['z_vmax'], 50)
        vmin = params['z_vmin']
        vmax = params['z_vmax']
        label = 'KuPR_c'
    elif choice==2:
        x = self.dpr.ds.zFactorMeasured[:,:,:,1].where(self.dpr.ds.zFactorMeasured[:,:,:,1] >= 15)
        cmap = 'Spectral_r'
        levels=np.linspace(params['z_vmin'], params['z_vmax'], 50)
        vmin = params['z_vmin']
        vmax = params['z_vmax']
        label = 'KaPR'
    elif choice==3:
        x = self.dpr.ds.zFactorFinal[:,:,:,1].where(self.dpr.ds.zFactorFinal[:,:,:,1] >= 15)
        cmap = 'Spectral_r'
        levels=np.linspace(params['z_vmin'], params['z_vmax'], 50)
        vmin = params['z_vmin']
        vmax = params['z_vmax']
        label = 'KaPR_c'
    elif choice==4:
        x = self.dpr.ds.zFactorMeasured[:,:,:,0].where(self.dpr.ds.zFactorMeasured[:,:,:,0] >= 10) - self.dpr.ds.zFactorMeasured[:,:,:,1].where(self.dpr.ds.zFactorMeasured[:,:,:,1] >= 15)
        cmap = cmaps.turbo
        levels=np.linspace(params['dfr_vmin'], params['dfr_vmax'], 50)
        vmin = params['dfr_vmin']
        vmax = params['dfr_vmax']
        label = 'DFR'
    elif choice==5:
        x = self.dpr.ds.zFactorFinal[:,:,:,0].where(self.dpr.ds.zFactorFinal[:,:,:,0] >= 10) - self.dpr.ds.zFactorFinal[:,:,:,1].where(self.dpr.ds.zFactorFinal[:,:,:,1] >= 15)
        cmap = cmaps.turbo
        levels=np.linspace(params['dfr_vmin'], params['dfr_vmax'], 50)
        vmin = params['dfr_vmin']
        vmax = params['dfr_vmax']
        label = 'DFR_c'
    elif choice==6:
        x = self.dpr.ds.paramDSD[:,:,:,1].where(self.dpr.ds.zFactorFinal[:,:,:,0] >= 10) 
        cmap = cmaps.magma
        levels=np.linspace(params['dm_vmin'], params['dm_vmax'], 50)
        vmin = params['dm_vmin']
        vmax = params['dm_vmax']
        label = 'Dm_2ADPR'
    elif choice==7:
        x = self.dpr.ds.paramDSD[:,:,:,0].where(self.dpr.ds.zFactorFinal[:,:,:,0] >= 10)/10 #convert to log(mm^-1 m^-3)
        cmap = cmaps.plasma
        levels=np.linspace(params['nw_vmin'], params['nw_vmax'], 50)
        vmin = params['nw_vmin']
        vmax = params['nw_vmax']
        label = 'Nw_2ADPR'
    elif choice==8:
        x = self.dpr.ds.precipRate.where(self.dpr.ds.zFactorFinal[:,:,:,0] >= 10) 
        x = x.where(x > 0)
        x = np.log10(x)
        cmap = cmaps.plasma
        levels=np.linspace(params['r_vmin'], params['r_vmax'], 50)
        vmin = params['r_vmin']
        vmax = params['r_vmax']
        label = 'R_2ADPR'
    elif choice==9:
        x = self.dpr.ds.Dml_nn
        cmap = cmaps.magma
        levels=np.linspace(params['dm_vmin'], params['dm_vmax'], 50)
        vmin = params['dm_vmin']
        vmax = params['dm_vmax']
        label = 'Dm_liq_NN'
    elif choice==10:
        x = self.dpr.ds.Dms_nn
        cmap = cmaps.magma
        levels=np.linspace(params['dm_vmin'], params['dm_vmax'], 50)
        vmin = params['dm_vmin']
        vmax = params['dm_vmax']
        label = 'Dm_sol_NN'
    elif choice==11:
        x = self.dpr.ds.R_nn
        x = x.where(x > 0)
        x = np.log10(x)
        cmap = cmaps.plasma
        levels=np.linspace(params['r_vmin'], params['r_vmax'], 50)
        vmin = params['r_vmin']
        vmax = params['r_vmax']
        label = 'R_NN'

    pm = ax.pcolormesh(self.dpr.ds.distance.values[s:e,w],self.dpr.ds.height.values[s,w,:]/1000,x.values[s:e,w,:].T,cmap=cmap,vmin=vmin,vmax=vmax,levels=levels)
    if params['temperature']:
        ax.contour(self.dpr.ds.distance.values[s:e,w],self.dpr.ds.height.values[s,w,:]/1000,self.dpr.ds.airTemperature.values[s:e,w,:].T-273.15,colors='k',levels=params['t_levels'])
    ax.set_ylim([0,params['y_max']])
    text = ax.text(0.025,0.85,label,fontsize=12,transform=ax.transAxes)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])   

def cross_section_method2(self,choice,along_track_index,ax,params): 
    from drpy.graph import cmaps
    if choice ==0:
        x = self.dpr.ds.zFactorMeasured[:,:,:,0].where(self.dpr.ds.zFactorMeasured[:,:,:,0] >= 10)
        cmap = 'Spectral_r'
        levels=np.linspace(params['z_vmin'], params['z_vmax'], 50)
        vmin = params['z_vmin']
        vmax = params['z_vmax']
        label = 'KuPR'
    elif choice==1:
        x = self.dpr.ds.zFactorFinal[:,:,:,0].where(self.dpr.ds.zFactorFinal[:,:,:,0] >= 10)
        cmap = 'Spectral_r'
        levels=np.linspace(params['z_vmin'], params['z_vmax'], 50)
        vmin = params['z_vmin']
        vmax = params['z_vmax']
        label = 'KuPR_c'
    elif choice==2:
        x = self.dpr.ds.zFactorMeasured[:,:,:,1].where(self.dpr.ds.zFactorMeasured[:,:,:,1] >= 15)
        cmap = 'Spectral_r'
        levels=np.linspace(params['z_vmin'], params['z_vmax'], 50)
        vmin = params['z_vmin']
        vmax = params['z_vmax']
        label = 'KaPR'
    elif choice==3:
        x = self.dpr.ds.zFactorFinal[:,:,:,1].where(self.dpr.ds.zFactorFinal[:,:,:,1] >= 15)
        cmap = 'Spectral_r'
        levels=np.linspace(params['z_vmin'], params['z_vmax'], 50)
        vmin = params['z_vmin']
        vmax = params['z_vmax']
        label = 'KaPR_c'
    elif choice==4:
        x = self.dpr.ds.zFactorMeasured[:,:,:,0].where(self.dpr.ds.zFactorMeasured[:,:,:,0] >= 10) - self.dpr.ds.zFactorMeasured[:,:,:,1].where(self.dpr.ds.zFactorMeasured[:,:,:,1] >= 15)
        cmap = cmaps.turbo
        levels=np.linspace(params['dfr_vmin'], params['dfr_vmax'], 50)
        vmin = params['dfr_vmin']
        vmax = params['dfr_vmax']
        label = 'DFR'
    elif choice==5:
        x = self.dpr.ds.zFactorFinal[:,:,:,0].where(self.dpr.ds.zFactorFinal[:,:,:,0] >= 10) - self.dpr.ds.zFactorFinal[:,:,:,1].where(self.dpr.ds.zFactorFinal[:,:,:,1] >= 15)
        cmap = cmaps.turbo
        levels=np.linspace(params['dfr_vmin'], params['dfr_vmax'], 50)
        vmin = params['dfr_vmin']
        vmax = params['dfr_vmax']
        label = 'DFR_c'
    elif choice==6:
        x = self.dpr.ds.paramDSD[:,:,:,1].where(self.dpr.ds.zFactorFinal[:,:,:,0] >= 10) 
        cmap = cmaps.magma
        levels=np.linspace(params['dm_vmin'], params['dm_vmax'], 50)
        vmin = params['dm_vmin']
        vmax = params['dm_vmax']
        label = 'Dm_2ADPR'
    elif choice==7:
        x = self.dpr.ds.paramDSD[:,:,:,0].where(self.dpr.ds.zFactorFinal[:,:,:,0] >= 10)/10 #convert to log(mm^-1 m^-3)
        cmap = cmaps.plasma
        levels=np.linspace(params['nw_vmin'], params['nw_vmax'], 50)
        vmin = params['nw_vmin']
        vmax = params['nw_vmax']
        label = 'Nw_2ADPR'
    elif choice==8:
        x = self.dpr.ds.precipRate.where(self.dpr.ds.zFactorFinal[:,:,:,0] >= 10) 
        x = x.where(x > 0)
        x = np.log10(x)
        cmap = cmaps.plasma
        levels=np.linspace(params['r_vmin'], params['r_vmax'], 50)
        vmin = params['r_vmin']
        vmax = params['r_vmax']
        label = 'R_2ADPR'
    elif choice==9:
        x = self.dpr.ds.Dml_nn
        cmap = cmaps.magma
        levels=np.linspace(params['dm_vmin'], params['dm_vmax'], 50)
        vmin = params['dm_vmin']
        vmax = params['dm_vmax']
        label = 'Dm_liq_NN'
    elif choice==10:
        x = self.dpr.ds.Dms_nn
        cmap = cmaps.magma
        levels=np.linspace(params['dm_vmin'], params['dm_vmax'], 50)
        vmin = params['dm_vmin']
        vmax = params['dm_vmax']
        label = 'Dm_sol_NN'
    elif choice==11:
        x = self.dpr.ds.R_nn
        x = x.where(x > 0)
        x = np.log10(x)
        cmap = cmaps.plasma
        levels=np.linspace(params['r_vmin'], params['r_vmax'], 50)
        vmin = params['r_vmin']
        vmax = params['r_vmax']
        label = 'R_NN'

    d = self.dpr.ds.distance.values[along_track_index,:]
    pm = ax.pcolormesh(np.tile(d[:,np.newaxis],(1,176)),self.dpr.ds.height.values[along_track_index,:,:]/1000,x.values[along_track_index,:,:],cmap=cmap,vmin=vmin,vmax=vmax,levels=levels)
    if params['temperature']:
        ax.contour(np.tile(d[:,np.newaxis],(1,176)),self.dpr.ds.height.values[along_track_index,:,:]/1000,self.dpr.ds.airTemperature.values[along_track_index,:,:]-273.15,colors='k',levels=params['t_levels'])
    ax.set_ylim([0,params['y_max']])
    text = ax.text(0.025,0.85,label,fontsize=12,transform=ax.transAxes)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])  