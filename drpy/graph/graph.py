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
        else:
            self.xrds = GPMDPR.xrds
        
        self.graphdict = {'title':None,'xlabel':None,'ylabel':None,'xlim':None,'ylim':None}
        
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