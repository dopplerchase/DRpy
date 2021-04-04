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

class APR_plot_obj():
    """Author: Randy J. Chase. This object has lots of built in plotting functions for the xarray dataset created by drpy.core.APR"""
    def __init__(self,APR=None): 
        if APR is None:
            self.xrds = None
        else:
            self.xrds = APR.xrds
            
    def threepanel(self,ylim=[0,10000],vlim=[-10,40],vlim2=[-2,10],scan=12):
    
        fig,axes = plt.subplots(3,1,figsize=(5,10))
        ax= axes[0]
        pm = ax.pcolormesh(self.xrds.time3d[:,scan,:],
                           self.xrds.alt3d[:,scan,:],self.xrds.Ku[:,scan,:],
                           cmap=cmaps.HomeyerRainbow,vmin=vlim[0],
                           vmax=vlim[1])
        
        cbar = plt.colorbar(pm,ax=ax)
        cbar.set_label('Z, [dBZ]')
        ax.set_ylabel('Alt, [m]')
        ax.set_title('Ku-band')
        ax.set_ylim([ylim[0],ylim[1]])
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))

        ax= axes[1]
        pm = ax.pcolormesh(self.xrds.time3d[:,scan,:],
                           self.xrds.alt3d[:,scan,:],self.xrds.Ka[:,scan,:],
                           cmap=cmaps.HomeyerRainbow,vmin=vlim[0],
                           vmax=vlim[1])
        
        cbar = plt.colorbar(pm,ax=ax)
        cbar.set_label('Z, [dBZ]')
        ax.set_ylabel('Alt, [m]')
        ax.set_title('Ka-band')
        ax.set_ylim([ylim[0],ylim[1]])
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))

        ax= axes[2]
        pm = ax.pcolormesh(self.xrds.time3d[:,scan,:],
                           self.xrds.alt3d[:,scan,:],self.xrds.Ku[:,scan,:]-self.xrds.Ka[:,scan,:],
                           cmap=cmaps.turbo,vmin=vlim2[0],
                           vmax=vlim2[1])
        
        cbar = plt.colorbar(pm,ax=ax)
        cbar.set_label('DFR, [dB]')
        ax.set_ylabel('Alt, [m]')
        ax.set_title('Ka-band')
        ax.set_ylim([ylim[0],ylim[1]])
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))

        plt.tight_layout()

        plt.show()


def make_colorbar(ax,vmin,vmax,cmap):
    cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap,
                              norm=cNorm,
                              orientation='horizontal',extend='both')
    return cb1

class case_study:

  def __init__(self,filename=None,center_lat=None,center_lon=None):

    import drpy 
    dpr = drpy.core.GPMDPR(filename=filename,outer_swath=True,auto_run=False)
    dpr.read()
    dpr.calc_heights()
    dpr.toxr(clutter=False,echotop=False)
    corners = [center_lon - 5,center_lon + 5,center_lat-5,center_lat+5]
    dpr.corners = corners
    dpr.setboxcoords()
    #drop dead weight (i.e. blank data)
    dpr.xrds = dpr.xrds.dropna(dim='along_track',how='all')
    self.dpr = dpr

  def plotter(self,start_index=25,end_index=-25,scan=24,params_new=None):

    if params_new is None:
      params = {'z_vmin':10,'z_vmax':40,'y_max':10}
    else:
      params = {'z_vmin':10,'z_vmax':40,'y_max':10}
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

    #determine center
    s = 0 
    e = self.dpr.xrds.lons.shape[0]
    middle = int((e-s)/2)
    lon0 = self.dpr.xrds.lons.values[middle,24]
    lat0 = self.dpr.xrds.lats.values[middle,24]
    #set specific aspect for cartopy plot
    x = 7.5
    y = 0.6666666*x
    corners = [lon0-x,lon0+x,lat0-y,lat0+y]

    #draw axes with proplot
    array = [  # the "picture" (1 == subplot A, 2 == subplot B, etc.)
    [1,1,1,0,0,0,0],
    [1,1,1, 4,4,4,4],
    [2,2,2,4,4,4,4],
    [2,2,2,4,4,4,4],
    [3,3,3,4,4,4,4],
    [3,3,3,0,0,0,0]]
    fig, axs = plot.subplots(array, width=10,height=5,span=False,proj={4:'cyl'},tight=False)
    fig.set_facecolor('w')

    #do things to map subplot 
    ax = axs[3]
    ax.format(borderscolor='k',gridminor=True,
              borders=True,lonlabels='t', latlabels='l',
              lonlim=(corners[0], corners[1]), latlim=(corners[2], corners[3]))
    
    box = ax.get_position()
    box.y0 = box.y0 - 0.05
    box.y1 = box.y1 - 0.05
    box.x0 = box.x0 + 0.01
    box.x1 = box.x1 + 0.01
    ax.set_position(box)
    # ax.set_extent(corners)
    ax.add_feature(cartopy.feature.LAND.with_scale('50m'),facecolor=[0.9,0.9,0.9])
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
    plt.setp(ax.spines.values(), color='orangered',lw=2)
    #plot swath (descending)
    ax.plot(self.dpr.xrds.lons[:,0]+0.0485,self.dpr.xrds.lats[:,0],'--k')
    ax.plot(self.dpr.xrds.lons[:,-1]-0.0485,self.dpr.xrds.lats[:,-1],'--k')

    self.dpr.extract_nearsurf()

    #plot data 
    pm = ax.scatter(self.dpr.xrds.lons[:,:],self.dpr.xrds.lats[:,:],c=self.dpr.xrds.nearsurfaceKu[:,:],vmin=params['z_vmin'],vmax=params['z_vmax'],s=1,cmap='Spectral_r',linewidths=0.1,zorder=10)
    s = start_index
    e = end_index
    w = scan

    ax.plot(self.dpr.xrds.lons[s:e,w],self.dpr.xrds.lats[s:e,w],'-k',markerfacecolor='w',ms=10,zorder=12)
    ax.plot(self.dpr.xrds.lons[s,w],self.dpr.xrds.lats[s,w],'k',markerfacecolor='w',ms=10,label='Start',marker='$L$',markeredgewidth=1,zorder=12)
    ax.plot(self.dpr.xrds.lons[e,w],self.dpr.xrds.lats[e,w],'k',markerfacecolor='w',ms=10,label='End',marker='$R$',markeredgewidth=1,zorder=12)

    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'),color='k',zorder=11)

    #add zoomed out map for context
    inset_axis = ax.inset([0.575,-0.4,0.5,0.5],proj='cyl',zoom=False,zorder=11)
    inset_axis.format(land=True,landcolor=[0.9,0.9,0.9],borderscolor='k',gridminor=True)
    inset_axis.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
    inset_axis.plot(self.dpr.xrds.lons[:,0]+0.0485,self.dpr.xrds.lats[:,0],'--k',lw=0.5,)
    inset_axis.plot(self.dpr.xrds.lons[:,-1]-0.0485,self.dpr.xrds.lats[:,-1],'--k',lw=0.5,)
    inset_axis.plot([corners[0],corners[0],corners[1],corners[1],corners[0]],[corners[2],corners[3],corners[3],corners[2],corners[2]],'-',color='orangered')
    timestr = pd.to_datetime(self.dpr.xrds.time[middle,24].values).strftime(format='%Y-%m-%d %H:%M')
    text = inset_axis.text(-0.05,-0.2,'Scan Time: ' + timestr,transform=ax.transAxes,fontsize=10)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
    text = inset_axis.text(0.025,-0.275,'Created with DRpy',transform=ax.transAxes,fontsize=10)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
    #draw colorbars in better spot
    ax_cbar = fig.add_axes([0.6, 0.93, 0.33, 0.015])
    text = ax_cbar.text(-0.25,0,'$Z_{e}$, [$dBZ$]',fontsize=10,transform=ax_cbar.transAxes)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
    cb1 = make_colorbar(ax_cbar,params['z_vmin'],params['z_vmax'],plt.cm.Spectral_r)
    ax_cbar = fig.add_axes([0.6, 0.85, 0.33, 0.015])
    text = ax_cbar.text(-0.33,0,'$DFR_{Ku-Ka}$, [$dB$]',fontsize=10,transform=ax_cbar.transAxes)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
    cb2 = make_colorbar(ax_cbar,-2,10,cmaps.turbo)

    #fill cross-sections 


    self.dpr.extract_nearsurf()
    self.dpr.get_physcial_distance(reference_point=[self.dpr.xrds.lons.values[s,w],self.dpr.xrds.lats.values[s,w]])

    ax = axs[0]
    ku = self.dpr.xrds.NSKu.where(self.dpr.xrds.NSKu >= 10)
    pm = ax.pcolormesh(self.dpr.xrds.distance.values[s:e,w],self.dpr.xrds.NSKu.alt.values[s,w,:],ku.values[s:e,w,:].T,cmap='Spectral_r',vmin=params['z_vmin'],vmax=params['z_vmax'],levels=np.linspace(params['z_vmin'], params['z_vmax'], 50))
    ax.set_ylim([0,params['y_max']])
    # ax.yaxis.set_ticks(np.arange(0,params['y_max']+2,2))
    ax.xaxis.set_ticklabels([])
    text = ax.text(0.025,0.85,'KuPR',fontsize=12,transform=ax.transAxes)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])

    ax = axs[1]
    ka = self.dpr.xrds.MSKa.where(self.dpr.xrds.MSKa >= 15)
    pm = ax.pcolormesh(self.dpr.xrds.distance.values[s:e,w],self.dpr.xrds.NSKu.alt.values[s,w,:],ka.values[s:e,w,:].T,cmap='Spectral_r',vmin=params['z_vmin'],vmax=params['z_vmax'],levels=np.linspace(params['z_vmin'], params['z_vmax'], 50))
    ax.set_ylim([0,params['y_max']])
    # ax.yaxis.set_ticks(np.arange(0,params['y_max']+2,2))
    ax.xaxis.set_ticklabels([])
    ax.set_ylabel('Alt ASL, [km]',labelpad=15)
    text = ax.text(0.025,0.85,'KaPR',fontsize=12,transform=ax.transAxes)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])

    ax = axs[2]
    pm = ax.pcolormesh(self.dpr.xrds.distance.values[s:e,w],self.dpr.xrds.MSKa.alt.values[s,w,:],ku.values[s:e,w,:].T-ka.values[s:e,w,:].T,cmap='turbo',vmin=-2,vmax=10,levels=np.linspace(-2, 10, 50))
    ax.set_ylim([0,params['y_max']])
    # ax.yaxis.set_ticks(np.arange(0,params['y_max']+2,2))
    ax.xaxis.set_ticklabels(ax.xaxis.get_ticklocs())
    text = ax.text(0.025,0.85,'DFR',fontsize=12,transform=ax.transAxes)
    text.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
    ax.set_xlabel('Along Track Distance, [km]',labelpad=8)

