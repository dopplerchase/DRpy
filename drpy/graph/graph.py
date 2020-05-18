import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
from . import cmaps

#plot parameters that I personally like, feel free to make these your own.
matplotlib.rcParams['axes.facecolor'] = [0.9,0.9,0.9]
matplotlib.rcParams['axes.labelsize'] = 14
matplotlib.rcParams['axes.titlesize'] = 14
matplotlib.rcParams['xtick.labelsize'] = 12
matplotlib.rcParams['ytick.labelsize'] = 12
matplotlib.rcParams['legend.fontsize'] = 12
matplotlib.rcParams['legend.facecolor'] = 'w'

class GPMDPR_plot_obj():
    """Author: Randy J. Chase. This object has lots of built in plotting functions for the xarray dataset created by drpy.core.GPMDPR"""

    def __init__(self,GPMDPR=None,notebook=False): 
        if GPMDPR is None:
            self.xrds = None
            self.corners=None
        else:
            self.xrds = GPMDPR.xrds
            self.corners = GPMDPR.corners
        
        self.graphdict = {'title':None,'xlabel':None,'ylabel':None,'xlim':None,'ylim':None}
        
    def mapper(self,z = None,extent=None,vmin=None,vmax=None,cmap='cividis',box=False):
        """
        This method creates a cartopy map of the data held in the xrds.
        z is what is desired to be plotted at each lat and lon pair. If 
        z is none, then it will default to just showing where all the footprints are. 
        
        """
        
        #import cartopy stuff
        import cartopy
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        import matplotlib.ticker as mticker
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
        import cartopy.io.shapereader as shpreader
        from cartopy.mpl.geoaxes import GeoAxes
        from mpl_toolkits.axes_grid1 import AxesGrid
        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
        
        #make figure
        fig = plt.figure(figsize=(10, 10))
        #add the map
        ax = fig.add_subplot(1, 1, 1,projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.STATES.with_scale('50m'),lw=0.5)
        ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
        ax.add_feature(cartopy.feature.LAND.with_scale('50m'), edgecolor='black',lw=0.5,facecolor=[0.95,0.95,0.95])
        ax.add_feature(cartopy.feature.LAKES.with_scale('50m'), edgecolor='black')
        ax.add_feature(cartopy.feature.RIVERS.with_scale('50m'))


        if (self.corners is not None) and box:
            ax.plot([corners[0],corners[0],corners[1],corners[1],corners[0]],[corners[2],corners[3],corners[3],corners[2],corners[2]],'-r',lw=3)
            ax.plot(center_lon,center_lat,'w*',ms=16,zorder=7)
            
        if z is None:
            ax.scatter(self.xrds.lons,self.xrds.lats,zorder=5)
        else:
            pm = ax.scatter(self.xrds.lons,self.xrds.lats,c=z,cmap=cmap,s=25,vmin=vmin,vmax=vmax,zorder=5)
            cbar = plt.colorbar(pm,ax=ax,shrink=0.5)
            cbar.set_label(z.name + ', [' + z.units + ']')
        
        ax.plot(self.xrds.lons,self.xrds.lats,'o',fillstyle='none',color='k',markeredgewidth=0.1,ms=4,zorder=6)
                
        if (self.corners is not None) and (extent is None):
            #for some reason set_extent crashes the session on colab. 
#             ax.set_extent(self.corners)
            ax.set_xlim([self.corners[0],self.corners[1]])
            ax.set_ylim([self.corners[2],self.corners[3]])
            ax.set_xticks(np.arange(self.corners[0], self.corners[1], 1), crs=ccrs.PlateCarree())
            ax.set_yticks(np.linspace(self.corners[2], self.corners[3], 5), crs=ccrs.PlateCarree())
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)
        elif (self.corners is not None) and (extent is not None):
            #for some reason set_extent crashes the session on colab.
#             ax.set_extent(extent)
            ax.set_xlim([extent[0],extent[1]])
            ax.set_ylim([extent[2],extent[3]])
            ax.set_xticks(np.linspace(extent[0], extent[1], 5), crs=ccrs.PlateCarree())
            ax.set_yticks(np.linspace(extent[2], extent[3], 5), crs=ccrs.PlateCarree())
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)
        elif (self.corners is None) and (extent is not None):
            #for some reason set_extent crashes the session on colab.
#             ax.set_extent(extent)
            ax.set_xlim([extent[0],extent[1]])
            ax.set_ylim([extent[2],extent[3]])
            ax.set_xticks(np.linspace(extent[0], extent[1], 5), crs=ccrs.PlateCarree())
            ax.set_yticks(np.linspace(extent[2], extent[3], 5), crs=ccrs.PlateCarree())
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)
            
        self.ax = ax 
            
    def CFAD(self,variablekey=None,bins=None,mincnt=20,graphdict=None):
        
        if self.xrds is None:
            print('ERROR, no GPMDPR object found. Please add one.')
        else:
            if graphdict is None:
                graphdict = self.graphdict

            if (bins is None and variablekey is None) or (bins is None and variablekey is 'NSKu'):
                #set some default ranges, these are for snow at Ku-band
                z_up = np.linspace(0,7,25)
                z_right = np.linspace(10,35,25)
                bins = [z_right,z_up]
                variablekey = 'NSKu'
                variable = self.xrds[variablekey].values.ravel()
                
            elif bins is None and variablekey is 'DFR':
                #set some default ranges, these are for snow at Ku-band
                z_up = np.linspace(0,7,25)
                z_right = np.linspace(-5,10,25)
                bins = [z_right,z_up]
                variable = self.xrds['NSKu'].values.ravel() - self.xrds['MSKa'].values.ravel()
            
            elif variablekey is not None:
                variable = self.xrds[variablekey].values.ravel()
                
            hist = np.histogram2d(variable,self.xrds.alt.values.ravel(),bins=bins)
            
            z_right = bins[0]
            z_up = bins[1]
            mat_tot = hist[0]
            mat_norm = np.zeros(mat_tot.shape)
            
            for i in np.arange(0,z_up.shape[0]-1):
                # row_sum = np.sum(mat_tot[:,i])
                row_sum = np.ma.max(mat_tot[:,i])
                mat_norm[:,i] = mat_tot[:,i]/row_sum


            midpoints_x = np.zeros(len(z_right)-1)
            for i in np.arange(0,len(z_right)-1):
                midpoints_x[i] = z_right[i] + (z_right[i+1] - z_right[i])/2.

            midpoints_y = np.zeros(len(z_up)-1)
            for i in np.arange(0,len(z_up)-1):
                midpoints_y[i] = z_up[i]  + (z_up[i+1] - z_up[i])/2.

            mat_norm= np.ma.masked_where(mat_tot < mincnt,mat_norm)
            fig = plt.figure(figsize=(5,5))
            fig.set_facecolor('w')
            ax = plt.gca()
            pm = ax.contourf(midpoints_x,midpoints_y,mat_norm.T,np.linspace(0.05,1,10),cmap=cmaps.turbo)
            ax.contour(midpoints_x,midpoints_y,mat_norm.T,np.linspace(0.05,1,10),colors='w',linewidths=0.2)


            cbar = plt.colorbar(pm,ax=ax)
            cbar.set_label('Normalized Frequency',fontsize=14)
            cbar.ax.tick_params(labelsize=12)
            ax.set_facecolor([0.9,0.9,0.9])
            ax.tick_params(labelsize=12)
            ax.set_xlabel(graphdict['xlabel'],fontsize=14)
            ax.set_ylabel(graphdict['ylabel'],fontsize=14)
            ax.set_title(graphdict['title'],fontsize=14)

            plt.tight_layout()
            
            
    def CFTD(self,variablekey=None,bins=None,mincnt=20,graphdict=None):
        
        if self.xrds is None:
            print('ERROR, no GPMDPR object found. Please add one.')
        else:

            if graphdict is None:
                graphdict = self.graphdict
                graphdict['ylabel'] = 'MERRA2 Temperature, [$\degree{C}$]'

            if bins is None and variablekey is None:
                #set some default ranges, these are for snow at Ku-band
                z_up = np.linspace(-25,5,25)
                z_right = np.linspace(10,35,25)
                bins = [z_right,z_up]
                variablekey = 'NSKu'
                variable = self.xrds[variablekey].values.ravel()
                #build label string 
                labelstr = variablekey + '  [' + self.xrds[variablekey].units + ']'
                graphdict['xlabel'] = labelstr
                
                
            elif bins is None and variablekey is 'DFR':
                #set some default ranges, these are for snow at Ku-band
                z_up = np.linspace(-25,5,25)
                z_right = np.linspace(-5,10,25)
                bins = [z_right,z_up]
                variable = self.xrds['NSKu'].values.ravel() - self.xrds['MSKa'].values.ravel()
                #build label string 
                labelstr = variablekey + ' [dB]'
                graphdict['xlabel'] = labelstr
                
            elif bins is None and variablekey is 'DFR_c':
                #set some default ranges, these are for snow at Ku-band
                z_up = np.linspace(-25,5,25)
                z_right = np.linspace(-5,10,25)
                bins = [z_right,z_up]
                variable = self.xrds['NSKu_c'].values.ravel() - self.xrds['MSKa_c'].values.ravel()
                #build label string 
                labelstr = variablekey + ' [dB]'
                graphdict['xlabel'] = labelstr
            
            elif variablekey is not None:
                variable = self.xrds[variablekey].values.ravel()
                #build label string 
                labelstr = variablekey + '  [' + self.xrds[variablekey].units + ']'
                graphdict['xlabel'] = labelstr
                
                if bins is None:
                    Q = np.nanpercentile(variable,[0.1,99.9])
                    z_up = np.linspace(-25,5,25)
                    z_right = np.linspace(Q[0],Q[1],25)
                    bins = [z_right,z_up]
                    
                    
                
            hist = np.histogram2d(variable,self.xrds.T.values.ravel()-273,bins=bins)
            
            z_right = bins[0]
            z_up = bins[1]
            mat_tot = hist[0]
            mat_norm = np.zeros(mat_tot.shape)
            
            for i in np.arange(0,z_up.shape[0]-1):
                # row_sum = np.sum(mat_tot[:,i])
                row_sum = np.ma.max(mat_tot[:,i])
                mat_norm[:,i] = mat_tot[:,i]/row_sum


            midpoints_x = np.zeros(len(z_right)-1)
            for i in np.arange(0,len(z_right)-1):
                midpoints_x[i] = z_right[i] + (z_right[i+1] - z_right[i])/2.

            midpoints_y = np.zeros(len(z_up)-1)
            for i in np.arange(0,len(z_up)-1):
                midpoints_y[i] = z_up[i]  + (z_up[i+1] - z_up[i])/2.

            mat_norm= np.ma.masked_where(mat_tot < mincnt,mat_norm)
            fig = plt.figure(figsize=(5,5))
            fig.set_facecolor('w')
            ax = plt.gca()
            pm = ax.contourf(midpoints_x,midpoints_y,mat_norm.T,np.linspace(0.05,1,10),cmap=cmaps.turbo)
            ax.contour(midpoints_x,midpoints_y,mat_norm.T,np.linspace(0.05,1,10),colors='w',linewidths=0.2)


            cbar = plt.colorbar(pm,ax=ax)
            cbar.set_label('Normalized Frequency',fontsize=14)
            cbar.ax.tick_params(labelsize=12)
            ax.set_facecolor([0.9,0.9,0.9])
            ax.tick_params(labelsize=12)
            ax.set_xlabel(graphdict['xlabel'],fontsize=14)
            ax.set_ylabel(graphdict['ylabel'],fontsize=14)
            ax.set_title(graphdict['title'],fontsize=14)
            ax.invert_yaxis()

            plt.tight_layout()
            
class APR_plot_obj():
    """Author: Randy J. Chase. This object has lots of built in plotting functions for the xarray dataset created by drpy.core.APR"""
    def __init__(self,APR=None): 
        if APR is None:
            self.xrds = None
        else:
            self.xrds = APR.xrds
            
    def threepanel(self,ylim=[0,10000],vlim=[-10,40],vlim2=[-2,10]):
    
        fig,axes = plt.subplots(3,1,figsize=(5,10))
        ax= axes[0]
        pm = ax.pcolormesh(self.xrds.time3d[:,12,:],
                           self.xrds.alt3d[:,12,:],self.xrds.Ku[:,12,:],
                           cmap=cmaps.HomeyerRainbow,vmin=vlim[0],
                           vmax=vlim[1])
        
        cbar = plt.colorbar(pm,ax=ax)
        cbar.set_label('Z, [dBZ]')
        ax.set_ylabel('Alt, [m]')
        ax.set_title('Ku-band')
        ax.set_ylim([ylim[0],ylim[1]])
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))

        ax= axes[1]
        pm = ax.pcolormesh(self.xrds.time3d[:,12,:],
                           self.xrds.alt3d[:,12,:],self.xrds.Ka[:,12,:],
                           cmap=cmaps.HomeyerRainbow,vmin=vlim[0],
                           vmax=vlim[1])
        
        cbar = plt.colorbar(pm,ax=ax)
        cbar.set_label('Z, [dBZ]')
        ax.set_ylabel('Alt, [m]')
        ax.set_title('Ka-band')
        ax.set_ylim([ylim[0],ylim[1]])
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))

        ax= axes[2]
        pm = ax.pcolormesh(self.xrds.time3d[:,12,:],
                           self.xrds.alt3d[:,12,:],self.xrds.Ku[:,12,:]-self.xrds.Ka[:,12,:],
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