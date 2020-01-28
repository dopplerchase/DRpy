from __future__ import absolute_import

import xarray as xr 
import h5py 
import numpy as np
import pandas as pd

class GPMDPR():
    """Author: Randy J. Chase. This object is intended to help with the efficient processing of GPM-DPR radar files. Currently, xarray cannot read the files directly. So here is an attempt to do so. Once in xarray format, the effcient search functions can be used. Currently, I do not have this function pass all variables through. Just the files I need to work with."""

    def __init__(self,filename=[],boundingbox=None): 
        self.filename = filename
        self.xrds = None
        self.datestr=None
        self.height= None
        self.corners = boundingbox
        self.retrieval_flag = 0
        self.interp_flag = 0
        
    def read(self):
        self.hdf = h5py.File(self.filename,'r')
        
        ###set some global parameters
        
        #whats the shape of the DPR files, mainly. 
        shape = self.hdf['NS']['PRE']['zFactorMeasured'][:,12:37,:].shape
        self.along_track = np.arange(0,shape[0])
        self.cross_track = np.arange(0,shape[1])
        self.range = np.arange(0,shape[2])
        
    def get_highest_clutter_bin(self):

        """
        This function makes us ground clutter conservative by supplying a clutter mask to apply to the fields 

        """

        ku = self.hdf['NS']['PRE']['binClutterFreeBottom'][:,12:37]
        ku = np.reshape(ku,[1,ku.shape[0],ku.shape[1]])
        ka = self.hdf['MS']['PRE']['binClutterFreeBottom'][:]
        ka = np.reshape(ka,[1,ka.shape[0],ka.shape[1]])
        both = np.vstack([ku,ka])
        pick_max = np.argmin(both,axis=0)
        ku = self.hdf['NS']['PRE']['binClutterFreeBottom'][:,12:37]
        ka = self.hdf['MS']['PRE']['binClutterFreeBottom'][:]
        inds_to_pick = np.zeros(ku.shape,dtype=int)
        ind = np.where(pick_max == 0)
        inds_to_pick[ind] = ku[ind]
        ind = np.where(pick_max == 1)
        inds_to_pick[ind] = ka[ind]

        dummy_matrix = np.ma.zeros([inds_to_pick.shape[0],inds_to_pick.shape[1],176])
        for i in np.arange(0,dummy_matrix.shape[0]):
            for j in np.arange(0,dummy_matrix.shape[1]):
                dummy_matrix[i,j,inds_to_pick[i,j]:] = 1

        self.dummy = np.ma.asarray(dummy_matrix,dtype=int)
        
    def echotop(self):
        keeper = self.range
        keeper = np.reshape(keeper,[1,keeper.shape[0]])
        keeper = np.tile(keeper,(25,1))
        keeper = np.reshape(keeper,[1,keeper.shape[0],keeper.shape[1]])
        keeper = np.tile(keeper,(self.xrds.NSKu_c.values.shape[0],1,1))
        keeper[np.isnan(self.xrds.NSKu_c)] = 9999

        inds_to_pick = np.argmin(keeper,axis=2)
        dummy_matrix = np.ma.zeros([inds_to_pick.shape[0],inds_to_pick.shape[1],176])
        for i in np.arange(0,dummy_matrix.shape[0]):
            for j in np.arange(0,dummy_matrix.shape[1]):
                dummy_matrix[i,j,:inds_to_pick[i,j]] = 1

        self.dummy2 = np.ma.asarray(dummy_matrix,dtype=int)
        
    def calcAltASL(self):
        
        """
        This function calculates the height of each radar gate above sea level. I am 
        not 100% this is exactly correct. Please use at your own risk!
        
        This is derived from some old code for TRMM (e.g., Stephen Nesbitt)
        
        """


        x2 = 2. * 17 #total degrees is 48 (from -17 to +17)
        re = 6378. #radius of the earth km 
        theta = -1 *(x2/2.) + (x2/48.)*np.arange(0,49) #break the -17 to 17 into equal degrees 

        theta2 = np.zeros(theta.shape[0]+1)
        theta = theta - 0.70833333/2. #shift thing to get left edge for pcolors
        theta2[:-1] = theta 
        theta2[-1] = theta[-1] + 0.70833333
        theta = theta2 * (np.pi/180.) #convert to radians

        prh = np.zeros([self.hdf['NS']['Longitude'][:,12:37].shape[0],49,176]) #set up matrix 
        for i in np.arange(0,176): #loop over num range gates
            for j in np.arange(0,49): #loop over scans 
                a = np.arcsin(((re+407)/re)*np.sin(theta[j]))-theta[j] #407 km is the orbit height, re radius of earth, 
                prh[:,j,i] = (176-(i))*0.125*np.cos(theta[j]+a) #more geometry 

        #reshape it into the same shape as the radar 
        self.height = prh[:,12:37,:]
        
    def toxr(self,snow=True,clutter=True,echotop=True):
    
        """
        
        This function converts the hdf file to an xarray dataset. 
        

        """
        #set the fillflag to false. 
        
        self.killflag = False
        #first thing first, check to make sure there are points in the bounding box.
        #cut points to make sure there are points in your box.This should save you time. 
        if self.corners is not None:
            
            #load data out of hdf
            lons = self.hdf['NS']['Longitude'][:,12:37]
            lats = self.hdf['NS']['Latitude'][:,12:37]
            #shove it into a dataarray
            da = xr.DataArray(self.hdf['MS']['Experimental']['flagSurfaceSnowfall'][:,:], dims=['along_track', 'cross_track'],
                           coords={'lons': (['along_track','cross_track'],lons),
                                   'lats': (['along_track','cross_track'],lats)})
            #cut the the edges of the box
            da = da.where((da.lons >= self.corners[0]) & (da.lons <= self.corners[1]) & (da.lats >= self.corners[2])  & (da.lats <= self.corners[3]),drop=False)
            #if you want only snowing cases, do the following
            if snow:
                da = da.where(da == 1,drop=False)
            #okay, now drop nans
            da = da.dropna(dim='along_track',how='all')
            #if there are no profiles, the len is 0, and we will set the kill flag
            if da.along_track.shape[0]==0:
                self.killflag = True
            
        #if there were no points it will not waste time with processing or io stuff    
        if self.killflag:
            pass
        else:          
            if self.datestr is None:
                self.parse_dtime()

            if self.height is None: 
                self.calcAltASL()

            da = xr.DataArray(self.hdf['MS']['Experimental']['flagSurfaceSnowfall'][:,:], dims=['along_track', 'cross_track'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr)})
            da.fillna(value=255)
            da.attrs['units'] = 'none'
            da.attrs['standard_name'] = 'experimental flag to diagnose snow at surface'

            #make xr dataset
            self.xrds = da.to_dataset(name = 'flagSurfaceSnow')
            #



            if clutter:
                self.get_highest_clutter_bin()
                da = xr.DataArray(self.dummy, dims=['along_track', 'cross_track','range'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr)})
                da.attrs['units'] = 'none'
                da.attrs['standard_name'] = 'flag to remove ground clutter'
                self.xrds['clutter'] = da

            da = xr.DataArray(self.hdf['NS']['SLV']['zFactorCorrectedNearSurface'][:,12:37], dims=['along_track', 'cross_track'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr)})
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'near surface Ku'
            da = da.where(da >= 12)
            self.xrds['nearsurfaceKu'] = da


            da = xr.DataArray(self.hdf['MS']['SLV']['zFactorCorrectedNearSurface'][:,:], dims=['along_track', 'cross_track'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats)})
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'near surface Ka'
            da = da.where(da >= 15)
            self.xrds['nearsurfaceKa'] = da

            da = xr.DataArray(self.hdf['NS']['SLV']['zFactorCorrected'][:,12:37,:], dims=['along_track', 'cross_track','range'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr),
                                               'alt':(['along_track', 'cross_track','range'],self.height)})
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'corrected KuPR'
            if clutter:
                da = da.where(self.xrds.clutter==0)
            da = da.where(da >= 12)
            self.xrds['NSKu_c'] = da

            da = xr.DataArray(self.hdf['MS']['SLV']['zFactorCorrected'][:,:,:], dims=['along_track', 'cross_track','range'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr),
                                               'alt':(['along_track', 'cross_track','range'],self.height)})
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'corrected KaPR, MS scan'
            if clutter:
                da = da.where(self.xrds.clutter==0)
            da = da.where(da >= 15)
            self.xrds['MSKa_c'] = da

            if echotop:
                self.echotop()
                da = xr.DataArray(self.dummy2, dims=['along_track', 'cross_track','range'],
                                           coords={'lons': (['along_track','cross_track'],lons),
                                                   'lats': (['along_track','cross_track'],lats),
                                                   'time': (['along_track','cross_track'],self.datestr)})
                da.attrs['units'] = 'none'
                da.attrs['standard_name'] = 'flag to remove noise outside cloud/precip top'
                self.xrds['echotop'] = da

            da = xr.DataArray(self.hdf['NS']['PRE']['zFactorMeasured'][:,12:37,:], dims=['along_track', 'cross_track','range'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr),
                                               'alt':(['along_track', 'cross_track','range'],self.height)})
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'measured KuPR'
            if clutter:
                da = da.where(self.xrds.clutter==0)
            if echotop:
                da = da.where(self.xrds.echotop==0)
            da = da.where(da >= 12)
            self.xrds['NSKu'] = da

            da = xr.DataArray(self.hdf['MS']['PRE']['zFactorMeasured'][:,:,:], dims=['along_track', 'cross_track','range'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr),
                                               'alt':(['along_track', 'cross_track','range'],self.height)})
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'measured KaPR, MS scan'
            if clutter:
                da = da.where(self.xrds.clutter==0)
            if echotop:
                da = da.where(self.xrds.echotop==0)
            da = da.where(da >= 15)
            self.xrds['MSKa'] = da



            da = xr.DataArray(self.hdf['NS']['SLV']['precipRate'][:,12:37,:], dims=['along_track', 'cross_track','range'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr),
                                               'alt':(['along_track', 'cross_track','range'],self.height)})
            da.attrs['units'] = 'mm hr^-1'
            da.attrs['standard_name'] = 'retrieved R, from DPR algo'
            if clutter:
                da = da.where(self.xrds.clutter==0)
            if echotop:
                da = da.where(self.xrds.echotop==0)
            self.xrds['R'] = da

            da = xr.DataArray(self.hdf['NS']['SLV']['paramDSD'][:,12:37,:,1], dims=['along_track', 'cross_track','range'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr),
                                               'alt':(['along_track', 'cross_track','range'],self.height)})
            da.attrs['units'] = 'mm'
            da.attrs['standard_name'] = 'retrieved Dm, from DPR algo'
            if clutter:
                da = da.where(self.xrds.clutter==0)
            if echotop:
                da = da.where(self.xrds.echotop==0)
            da = da.where(da >= 0)
            self.xrds['Dm_dpr'] = da


            if snow:
                self.xrds = self.xrds.where(self.xrds.flagSurfaceSnow==1,drop=False)
                
            if self.corners is not None:
                self.setboxcoords()
            
            
         
    def setboxcoords(self):
        
        if len(self.corners) > 0:
            self.ll_lon = self.corners[0]
            self.ur_lon = self.corners[1]
            self.ll_lat = self.corners[2]
            self.ur_lat = self.corners[3]
            self.xrds = self.xrds.where((self.xrds.lons >= self.ll_lon) & (self.xrds.lons <= self.ur_lon) & (self.xrds.lats >= self.ll_lat)  & (self.xrds.lats <= self.ur_lat),drop=False)
        else:
            pass
        
    def parse_dtime(self):
        
        """
        
        create datetime objects from the hdf file, typically run this after you already filtered for precip/snow
        to save computation time. 
        
        """
        year = np.asarray(self.hdf['MS']['ScanTime']['Year'][:],dtype=str)
        year = list(year)

        month = self.hdf['MS']['ScanTime']['Month'][:]
        ind = np.where(month < 10)
        month = np.asarray(month,dtype=str)
        month = numpy.char.rjust(month, 2, fillchar='0')
        month = list(month)

        day = self.hdf['MS']['ScanTime']['DayOfMonth'][:]
        ind = np.where(day < 10)
        day = np.asarray(day,dtype=str)
        day = numpy.char.rjust(day, 2, fillchar='0')
        day = list(day)

        hour = self.hdf['MS']['ScanTime']['Hour'][:]
        ind = np.where(hour < 10)
        hour = np.asarray(hour,dtype=str)
        hour = numpy.char.rjust(hour, 2, fillchar='0')
        hour = list(hour)

        minute = self.hdf['MS']['ScanTime']['Minute'][:]
        ind = np.where(minute < 10)
        minute = np.asarray(minute,dtype=str)
        minute = numpy.char.rjust(minute, 2, fillchar='0')
        minute = list(minute)

        second = self.hdf['MS']['ScanTime']['Second'][:]
        ind = np.where(second < 10)
        second = np.asarray(second,dtype=str)
        second = numpy.char.rjust(second, 2, fillchar='0')
        second = list(second)

        datestr  = [year[i] +"-"+ month[i]+ "-" + day[i] + ' ' + hour[i] + ':' + minute[i] + ':' + second[i]  for i in range(len(year))]
        datestr = np.asarray(datestr,dtype=str)
        datestr = np.reshape(datestr,[len(datestr),1])
        datestr = np.tile(datestr,(1,25))
        
        self.datestr = np.asarray(datestr,dtype=np.datetime64)
    
    def run_retrieval(self,path_to_models=None):
        
        if path_to_models is None:
            print('Please insert path to NN models')
        else:
            scaler = joblib.load(path_to_models+'scaler_training_unrimed.save')
            model_dual_dm = tf.keras.models.load_model(path_to_models+'NN_ZKuKa_Dm_deep_unrimed.h5',
                            custom_objects=None,
                            compile=True)
            
            #now we have to reshape things to make sure they are in the right shape for the NN model [n_samples,n_features]
            Ku = self.xrds.NSKu.values
            shape_step1 = Ku.shape
            Ku = Ku.reshape([Ku.shape[0],Ku.shape[1]*Ku.shape[2]])
            shape_step2 = Ku.shape
            Ku = Ku.reshape([Ku.shape[0]*Ku.shape[1]])
            Ka = self.xrds.MSKa.values
            Ka = Ka.reshape([Ka.shape[0],Ka.shape[1]*Ka.shape[2]])
            Ka = Ka.reshape([Ka.shape[0]*Ka.shape[1]])
            
            #Make sure we only run in on non-nan values. 
            ind_masked = np.isnan(Ku) 
            ind_masked2 = np.isnan(Ka)
            Ku_nomask = np.zeros(Ku.shape)
            Ka_nomask = np.zeros(Ka.shape)
            Ku_nomask[~ind_masked] = Ku[~ind_masked]
            Ka_nomask[~ind_masked] = Ka[~ind_masked]
            
            ind = np.where(Ku_nomask!=0)[0]

            #scale the input vectors by the mean that it was trained with
            X = np.zeros([Ku_nomask.shape[0],2])
            X[:,0] = (10**(Ku_nomask/10) - scaler.mean_[0])/scaler.scale_[0]
            X[:,1] = (10**(Ka_nomask/10)- scaler.mean_[2])/scaler.scale_[2]
            #
            
            #actually run the retrieval
            Dm_retrieved = model_dual_dm.predict(X[ind,:],batch_size=len(X[ind,0])) #dual
            
            #shove it back into the original shape
            Dm = np.zeros(Ku_nomask.shape)
            Dm[ind] = np.squeeze(Dm_retrieved)
            Dm = Dm.reshape([shape_step2[0],shape_step2[1]])
            Dm = Dm.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
            
            da = xr.DataArray(Dm, dims=['along_track', 'cross_track','range'],
                           coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                                   'lats': (['along_track','cross_track'],self.xrds.lons),
                                   'time': (['along_track','cross_track'],self.xrds.time),
                                   'alt':(['along_track', 'cross_track','range'],self.xrds.alt)})
            da.attrs['units'] = 'mm'
            da.attrs['standard_name'] = 'retrieved Dm from the NN (Chase et al. 2020)'
            da = da.where(da > 0.)
            self.xrds['Dm'] = da
            
            self.retrieval_flag = 1
        
    def get_merra(self,interp=True):
        
        """ 
        
        This is specific to my file structure at UIUC. Buyer beware.
        
        Future hope: interpolate merra-2 to the GPM-DPR vertical bins 
        
        """
        
        time = self.xrds.time.values
        orig_shape = time.shape
        time = np.reshape(time,[orig_shape[0]*orig_shape[1]])
        dates = pd.to_datetime(time,infer_datetime_format=True)
        dates = dates.to_pydatetime()
        dates = np.reshape(dates,[orig_shape[0],orig_shape[1]])

        year = dates[0,0].year
        month = dates[0,0].month
        day = dates[0,0].day

        if month < 10:
            month = '0'+ str(month)
        else:
            month = str(month)

        if day <10:
            day = '0' + str(day)
        else:
            day = str(day)

        ds_url = '/data/gpm/a/randyjc2/MERRA/NEW/'+ str(year) + '/' + 'MERRA2_400.inst6_3d_ana_Np.'+ str(year) + month + day+ '.nc4'

        ###load file
        merra = xr.open_dataset(ds_url)
        ###

        #select the closest profile to the lat, lon, time
        sounding = merra.sel(lon=self.xrds.lons,lat=self.xrds.lats,time=self.xrds.time,method='nearest')

        
        self.sounding = sounding
        
        if interp:
            self.interp_MERRA(keyname='T')
            self.interp_MERRA(keyname='U')
            self.interp_MERRA(keyname='V')
            self.interp_MERRA(keyname='QV')
            self.interp_flag = 1
            
        merra.close()
        
    def interp_MERRA(self,keyname=None):
        
        """ 
        
        This interpolates the MERRA two variables to the same veritcal levels as the GPM-DPR
        
        NOTE: I am not sure this is optimized! Not very fast..., but if you want you can turn it off
        
        """

        H_Merra = self.sounding.H.values
        H_gpm = self.xrds.alt.values
        new_variable = np.zeros(H_gpm.shape)
        for i in self.sounding.along_track.values:
            for j in self.sounding.cross_track.values:
                #fit func
                da = xr.DataArray(self.sounding[keyname].values[i,j,:], [('height', H_Merra[i,j,:]/1000)])
                da = da.interp(height=H_gpm[i,j,:])
                new_variable[i,j,:] = da.values


        da = xr.DataArray(new_variable, dims=['along_track', 'cross_track','range'],
               coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                       'lats': (['along_track','cross_track'],self.xrds.lons),
                       'time': (['along_track','cross_track'],self.xrds.time),
                       'alt':(['along_track', 'cross_track','range'],self.xrds.alt)})

        da.attrs['units'] = self.sounding[keyname].units
        da.attrs['standard_name'] = 'Interpolated ' + self.sounding[keyname].standard_name + ' to GPM height coord'
        self.xrds[keyname] = da
        
    def extract_nearsurf(self):
        keeper = self.xrds.range.values
        keeper = np.reshape(keeper,[1,keeper.shape[0]])
        keeper = np.tile(keeper,(25,1))
        keeper = np.reshape(keeper,[1,keeper.shape[0],keeper.shape[1]])
        keeper = np.tile(keeper,(self.xrds.NSKu.values.shape[0],1,1))
        keeper[np.isnan(self.xrds.NSKu.values)] = -9999

        inds_to_pick = np.argmax(keeper,axis=2)
        dummy_matrix = np.ma.zeros([inds_to_pick.shape[0],inds_to_pick.shape[1],176])

        #note, for all nan columns, it will say its 0, or the top of the GPM index, which should alway be nan anyway
        for i in np.arange(0,dummy_matrix.shape[0]):
            for j in np.arange(0,dummy_matrix.shape[1]):
                dummy_matrix[i,j,inds_to_pick[i,j]] = 1

        self.lowest_gate_index = np.ma.asarray(dummy_matrix,dtype=int)

        self.grab_variable(keyname='NSKu')
        self.grab_variable(keyname='NSKu_c')
        self.grab_variable(keyname='MSKa')
        self.grab_variable(keyname='MSKa_c')
        self.grab_variable(keyname='R')
        self.grab_variable(keyname='Dm_dpr')
        self.grab_variable(keyname='alt')
        
        if self.retrieval_flag == 1:
            self.grab_variable(keyname='Dm')

        if self.interp_flag == 1:
            self.grab_variable(keyname='T')
            self.grab_variable(keyname='U')
            self.grab_variable(keyname='V')
            self.grab_variable(keyname='QV')
        

    def grab_variable(self,keyname=None):

        if keyname is None:
            print('please supply keyname')
        else:
            variable = np.zeros([self.xrds.along_track.shape[0],self.xrds.cross_track.values.shape[0]])
            ind = np.where(self.lowest_gate_index == 0)
            variable[ind[0],ind[1]] = np.nan
            ind = np.where(self.lowest_gate_index == 1)
            variable[ind[0],ind[1]] = self.xrds[keyname].values[ind[0],ind[1],ind[2]]
            da = xr.DataArray(variable, dims=['along_track', 'cross_track'],
                       coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                               'lats': (['along_track','cross_track'],self.xrds.lons),
                               'time': (['along_track','cross_track'],self.xrds.time)})
            da = da.where(self.xrds.flagSurfaceSnow == 1)
            if keyname=='alt':
                da.attrs['units'] = 'km'
                da.attrs['standard_name'] = 'altitude of the near-surface bin'
                self.xrds[keyname+'_nearSurf'] = da
            else:
                da.attrs['units'] = self.xrds[keyname].units
                da.attrs['standard_name'] = 'near-surface' + self.xrds[keyname].standard_name
                self.xrds[keyname+'_nearSurf'] = da