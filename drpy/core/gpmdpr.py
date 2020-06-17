from __future__ import absolute_import
import xarray as xr 
import h5py 
import numpy as np
import pandas as pd
import datetime
import scipy
import scipy.interpolate

#turn off warnings so i can use the progressbar
import warnings
warnings.filterwarnings('ignore')



class GPMDPR():
    """
    Author: Randy J. Chase. This class is intended to help with the efficient processing of GPM-DPR radar files. 
    Currently, xarray cannot read NASA's HDF files directly (2A.GPM.DPR*). So here is an attempt to do so. 
    Once in xarray format, the effcient search functions can be used. 
    
    **NOTE: Currently, I do not have this function pass all variables through (there is quite the list of them.
    Maybe in the future I will generalize it to do so. But right now its a bit tedious to code up all the units and such
    
    Feel free to reach out to me on twitter (@dopplerchase) or email randyjc2@illinois.edu
    
    For your reference, please check out the ATBD: https://pps.gsfc.nasa.gov/GPMprelimdocs.html 
    """

    def __init__(self,filename=[],boundingbox=None): 
        self.filename = filename
        self.xrds = None
        self.datestr=None
        self.height= None
        self.corners = boundingbox
        self.retrieval_flag = 0
        self.interp_flag = 0
        
    def read(self):
        """
        This method simply reads the HDF file and gives it to the class. 
        
        """
        
        self.hdf = h5py.File(self.filename,'r')
        
        ###set some global parameters
        #whats the common shape of the DPR files
        shape = self.hdf['NS']['PRE']['zFactorMeasured'][:,12:37,:].shape
        self.along_track = np.arange(0,shape[0])
        self.cross_track = np.arange(0,shape[1])
        self.range = np.arange(0,shape[2])
        
    def get_highest_clutter_bin(self):
        """
        This method makes us ground clutter conservative by supplying a clutter mask to apply to the fields.
        It is based off the algorithim output of 'binClutterFreeBottom', which can be a bit conservative (~ 1km)

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
        """
        This method takes the already clutter filtered data for the corrected reflectivity and cuts the 
        noisy uncorrected reflectivity to the same height. Again, the method is a bit conservative, but is 
        a good place to start.
        
        """
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
        
    def toxr(self,ptype=None,clutter=True,echotop=True):
        """
        This is the main method of the package. It directly creates the xarray dataset from the HDF file. 
        
        To save computational time, it does first check to see if you set a box of interest. 
        Then it uses xarray effcient searching to make sure there are some profiles in that box. 
        

        """
        #set the precip type of interest. If none, give back all data...
        self.ptype= ptype
        self.snow = False
        self.precip = False
        

        if (self.ptype=='precip') or (self.ptype=='Precip') or (self.ptype=='PRECIP') or (self.ptype=='snow') or (self.ptype=='Snow') or (self.ptype=='SNOW'):
            self.precip=True
            if (self.ptype=='snow') or (self.ptype=='Snow') or (self.ptype=='SNOW'):
                self.snow=True
        
        #set the killflag to false. If this is True at the end, it means no points in the box were found. 
        self.killflag = False
        #first thing first, check to make sure there are points in the bounding box.
        #cut points to make sure there are points in your box.This should save you time. 
        if self.corners is not None:
            #load data out of hdf
            lons = self.hdf['NS']['Longitude'][:,12:37]
            lats = self.hdf['NS']['Latitude'][:,12:37]
            #shove it into a dataarray
            da = xr.DataArray(np.zeros(lons.shape), dims=['along_track', 'cross_track'],
                           coords={'lons': (['along_track','cross_track'],lons),
                                   'lats': (['along_track','cross_track'],lats)})
            #cut the the edges of the box
            da = da.where((da.lons >= self.corners[0]) & (da.lons <= self.corners[1]) & (da.lats >= self.corners[2])  & (da.lats <= self.corners[3]),drop=False)
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
                height = xr.open_dataarray('/data/gpm/a/randyjc2/DRpy/drpy/models/HEIGHTS.nc')
                height = height.values[np.newaxis,:,:]
                height = np.tile(height,(self.hdf['NS']['Longitude'].shape[0],1,1))
                self.height = height
                
            if self.corners is None:
                #load data out of hdf
                lons = self.hdf['NS']['Longitude'][:,12:37]
                lats = self.hdf['NS']['Latitude'][:,12:37]
            
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
            
            #ADD BBtop and Bottom 
            da = xr.DataArray(self.hdf['NS']['CSF']['binBBTop'][:,12:37], dims=['along_track', 'cross_track'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr)})
            da.fillna(value=-9999)
            da.attrs['units'] = 'none'
            da.attrs['standard_name'] = 'ind of BBtop'
            self.xrds['binBBTop'] = da
            
            da = xr.DataArray(self.hdf['NS']['CSF']['binBBBottom'][:,12:37], dims=['along_track', 'cross_track'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr)})
            da.fillna(value=-9999)
            da.attrs['units'] = 'none'
            da.attrs['standard_name'] = 'ind of BBtop'
            self.xrds['binBBBottom'] = da
            

            da = xr.DataArray(self.hdf['MS']['PRE']['flagPrecip'][:,:], dims=['along_track', 'cross_track'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr)})
            da.fillna(value=-9999)
            da.attrs['units'] = 'none'
            da.attrs['standard_name'] = 'flag to diagnose precip at surface. 11 is precip from both, 10 is preicp from just Ku-band'

            #make xr dataset
            self.xrds['flagPrecip'] = da
            #
            
            typePrecip = self.hdf['MS']['CSF']['typePrecip'][:]
            typePrecip = np.asarray(typePrecip,dtype=float)
            ind = np.where(typePrecip == -1111)
            typePrecip[ind] = np.nan
            ind = np.where(typePrecip == -9999)
            typePrecip[ind] = np.nan

            typePrecip = np.trunc(typePrecip/10000000)
            typePrecip = np.asarray(typePrecip,dtype=int)

            da = xr.DataArray(typePrecip, dims=['along_track', 'cross_track'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr)})
            da.fillna(value=-9999)
            da.attrs['units'] = 'none'
            da.attrs['standard_name'] = 'flag to diagnose raintype. If 1: Strat. If 2: Conv. If 3:other '
            
            self.xrds['typePrecip'] = da


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
            
            da = xr.DataArray(self.hdf['NS']['SLV']['epsilon'][:,12:37,:],dims=['along_track', 'cross_track','range'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr),
                                               'alt':(['along_track', 'cross_track','range'],self.height)})
            da.fillna(value=-9999.9)
            da = da.where(da >= 0)
            da.attrs['units'] = 'none'
            da.attrs['standard_name'] = 'epsilon value for retrieval'
            self.xrds['epsilon'] = da

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
            
            
            da = xr.DataArray(self.hdf['NS']['SLV']['paramDSD'][:,12:37,:,0], dims=['along_track', 'cross_track','range'],
                                       coords={'lons': (['along_track','cross_track'],lons),
                                               'lats': (['along_track','cross_track'],lats),
                                               'time': (['along_track','cross_track'],self.datestr),
                                               'alt':(['along_track', 'cross_track','range'],self.height)})
            da.attrs['units'] = 'dBNw'
            da.attrs['standard_name'] = 'retrieved Nw, from DPR algo'
            if clutter:
                da = da.where(self.xrds.clutter==0)
            if echotop:
                da = da.where(self.xrds.echotop==0)
            da = da.where(da >= 0)
            self.xrds['Nw_dpr'] = da

            if self.precip:
                #change this to 10 if you want to relax the conditions, because the ka band has bad sensativity
                self.xrds = self.xrds.where(self.xrds.flagPrecip==11)
#                 if self.snow:
#                     self.xrds = self.xrds.where(self.xrds.flagSurfaceSnow==1)
                    
            if self.corners is not None:
                self.setboxcoords()
            
            #to reduce size of data, drop empty cross-track sections 
#             self.xrds = self.xrds.dropna(dim='along_track',how='all')
            
            #as before, makes sure there is data...
            if self.xrds.along_track.shape[0]==0:
                self.killflag = True
            
            
         
    def setboxcoords(self):
        """
        This method sets all points outside the box to nan. 
        """
        
        if len(self.corners) > 0:
            self.ll_lon = self.corners[0]
            self.ur_lon = self.corners[1]
            self.ll_lat = self.corners[2]
            self.ur_lat = self.corners[3]
            self.xrds = self.xrds.where((self.xrds.lons >= self.ll_lon) & (self.xrds.lons <= self.ur_lon) & (self.xrds.lats >= self.ll_lat)  & (self.xrds.lats <= self.ur_lat),drop=False)
        else:
            print('ERROR, not boxcoods set...did you mean to do this?')
        
    def parse_dtime(self):
        """
        This method creates datetime objects from the hdf file in a timely mannor.
        Typically run this after you already filtered for precip/snow to save additional time. 
        """
        year = self.hdf['MS']['ScanTime']['Year'][:]
        ind = np.where(year == -9999)[0]
        year = np.asarray(year,dtype=str)
        year = list(year)

        month = self.hdf['MS']['ScanTime']['Month'][:]
        month = np.asarray(month,dtype=str)
        month = np.char.rjust(month, 2, fillchar='0')
        month = list(month)

        day = self.hdf['MS']['ScanTime']['DayOfMonth'][:]
        day = np.asarray(day,dtype=str)
        day = np.char.rjust(day, 2, fillchar='0')
        day = list(day)

        hour = self.hdf['MS']['ScanTime']['Hour'][:]
        hour = np.asarray(hour,dtype=str)
        hour = np.char.rjust(hour, 2, fillchar='0')
        hour = list(hour)

        minute = self.hdf['MS']['ScanTime']['Minute'][:]
        minute = np.asarray(minute,dtype=str)
        minute = np.char.rjust(minute, 2, fillchar='0')
        minute = list(minute)

        second = self.hdf['MS']['ScanTime']['Second'][:]
        second = np.asarray(second,dtype=str)
        second = np.char.rjust(second, 2, fillchar='0')
        second = list(second)
        
        datestr  = [year[i] +"-"+ month[i]+ "-" + day[i] + ' ' + hour[i] + ':' + minute[i] + ':' + second[i]  for i in range(len(year))]
        datestr = np.asarray(datestr,dtype=str)
        datestr[ind] = '1970-01-01 00:00:00'
        datestr = np.reshape(datestr,[len(datestr),1])
        datestr = np.tile(datestr,(1,25))

        self.datestr = np.asarray(datestr,dtype=np.datetime64)
    
    def run_retrieval(self,path_to_models=None,old=False,notebook=False):
        
        """
        This method is a way to run our neural network trained retreival to get Dm in snowfall. 
        Please see this AMS presentation until the paper comes out: *LINK HERE*.
        
        This method requires the use of tensorflow. So go install that. 
        
        """
        #load scalers
        from pickle import load
        if old:
            scaler_X = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/scaler_X.pkl', 'rb'))
            scaler_y = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/scaler_y.pkl', 'rb'))
        else:
            scaler_X = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/scaler_X_V2.pkl', 'rb'))
            scaler_y = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/scaler_y_V2.pkl', 'rb'))

        import tensorflow as tf
        from tensorflow.python.keras import losses
        
        #set number of threads = 1, this was crashing my parallel code
        if notebook:
            pass
        else:
            tf.config.threading.set_inter_op_parallelism_threads(1)
#         tf.config.threading.set_intra_op_parallelism_threads(1)
#         print('Number of threads set to {}'.format(tf.config.threading.get_inter_op_parallelism_threads()))


        

        
        #supress warnings. skrews up my progress bar when running in parallel
        def warn(*args, **kwargs):
            pass
        import warnings
        warnings.warn = warn
        

        
        if path_to_models is None:
            print('Please insert path to NN models')
        else:
            if old:
                model = tf.keras.models.load_model(path_to_models + 'NN_4by8.h5',custom_objects=None,compile=True)
            else:
                model = tf.keras.models.load_model(path_to_models + 'NN_6by8.h5',custom_objects=None,compile=True)
                
            #now we have to reshape things to make sure they are in the right shape for the NN model [n_samples,n_features]
            Ku = self.xrds.NSKu.values
            shape_step1 = Ku.shape
            Ku = Ku.reshape([Ku.shape[0],Ku.shape[1]*Ku.shape[2]])
            shape_step2 = Ku.shape
            Ku = Ku.reshape([Ku.shape[0]*Ku.shape[1]])
            Ka = self.xrds.MSKa.values
            Ka = Ka.reshape([Ka.shape[0],Ka.shape[1]*Ka.shape[2]])
            Ka = Ka.reshape([Ka.shape[0]*Ka.shape[1]])
            
            T = self.xrds['T'].values - 273.15 #expects in degC
            T = T.reshape([T.shape[0],T.shape[1]*T.shape[2]])
            T = T.reshape([T.shape[0]*T.shape[1]])
            
            #Make sure we only run in on non-nan values. 
            ind_masked = np.isnan(Ku) 
            ind_masked2 = np.isnan(Ka)
            Ku_nomask = np.zeros(Ku.shape)
            Ka_nomask = np.zeros(Ka.shape)
            T_nomask = np.zeros(T.shape)
            Ku_nomask[~ind_masked] = Ku[~ind_masked]
            Ka_nomask[~ind_masked] = Ka[~ind_masked]
            T_nomask[~ind_masked] = T[~ind_masked]
            
            ind = np.where(Ku_nomask!=0)[0]

            #scale the input vectors by the mean that it was trained with
            X = np.zeros([Ku_nomask.shape[0],3])
            X[:,0] = (Ku_nomask - scaler_X.mean_[0])/scaler_X.scale_[0] #ku 
            X[:,1] = ((Ku_nomask - Ka_nomask)- scaler_X.mean_[1])/scaler_X.scale_[1] #dfr
            X[:,2] = (T_nomask - scaler_X.mean_[2])/scaler_X.scale_[2] #T
            #
            
            yhat = model.predict(X[ind,0:3],batch_size=len(X[ind,0]))
            yhat = scaler_y.inverse_transform(yhat)
            yhat[:,1] = 10**yhat[:,1] #unlog Dm liquid
            yhat[:,2] = 10**yhat[:,2] #unlog Dm solid
            
            ind = np.where(Ku_nomask!=0)[0]
            Nw = np.zeros(Ku_nomask.shape)
            Nw[ind] = np.squeeze(yhat[:,0])
            Nw = Nw.reshape([shape_step2[0],shape_step2[1]])
            Nw = Nw.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
            Nw = np.ma.masked_where(Nw==0.0,Nw)

            da = xr.DataArray(Nw, dims=['along_track', 'cross_track','range'],
               coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                       'lats': (['along_track','cross_track'],self.xrds.lons),
                       'time': (['along_track','cross_track'],self.xrds.time),
                       'alt':(['along_track', 'cross_track','range'],self.xrds.alt)})
            
            da.fillna(value=-9999)
            da.attrs['units'] = 'log(m^-4)'
            da.attrs['standard_name'] = 'retrieved Nw from the NN (Chase et al. 2020)'
            da = da.where(da > 0.)
            self.xrds['Nw'] = da
           
            Dm = np.zeros(Ku_nomask.shape)
            Dm[ind] = np.squeeze(yhat[:,1])
            Dm = Dm.reshape([shape_step2[0],shape_step2[1]])
            Dm = Dm.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
            Dm = np.ma.masked_where(Dm==0.0,Dm)

            da = xr.DataArray(Dm, dims=['along_track', 'cross_track','range'],
               coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                       'lats': (['along_track','cross_track'],self.xrds.lons),
                       'time': (['along_track','cross_track'],self.xrds.time),
                       'alt':(['along_track', 'cross_track','range'],self.xrds.alt)})

            da.fillna(value=-9999)
            da.attrs['units'] = 'mm'
            da.attrs['standard_name'] = 'retrieved Liquid Eq. Dm from the NN (Chase et al. 2020)'

            self.xrds['Dm'] = da

            Dm_frozen = np.zeros(Ku_nomask.shape)
            Dm_frozen[ind] = np.squeeze(yhat[:,2])
            Dm_frozen = Dm_frozen.reshape([shape_step2[0],shape_step2[1]])
            Dm_frozen = Dm_frozen.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
            Dm_frozen = np.ma.masked_where(Dm_frozen==0.0,Dm_frozen)

            da = xr.DataArray(Dm_frozen, dims=['along_track', 'cross_track','range'],
               coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                       'lats': (['along_track','cross_track'],self.xrds.lons),
                       'time': (['along_track','cross_track'],self.xrds.time),
                       'alt':(['along_track', 'cross_track','range'],self.xrds.alt)})

            da.fillna(value=-9999)
            da.attrs['units'] = 'mm'
            da.attrs['standard_name'] = 'retrieved Frozen Dm from the NN (Chase et al. 2020)'

            self.xrds['Dm_frozen'] = da
        
            Nw = 10**Nw #undo log, should be in m^-4
            Dm = Dm/1000. # convert to m ^4
            IWC = (Nw*(Dm)**4*1000*np.pi)/4**(4) # the 1000 is density of water (kg/m^3)
            IWC = IWC*1000 #convert to g/m^3 
        
            da = xr.DataArray(IWC, dims=['along_track', 'cross_track','range'],
               coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                       'lats': (['along_track','cross_track'],self.xrds.lons),
                       'time': (['along_track','cross_track'],self.xrds.time),
                       'alt':(['along_track', 'cross_track','range'],self.xrds.alt)})

            da.fillna(value=-9999)
            da.attrs['units'] = 'g m^{-3}'
            da.attrs['standard_name'] = 'Calc IWC from retrieved Nw and Dm from the NN (Chase et al. 2020)'

            self.xrds['IWC'] = da
            
            self.retrieval_flag = 1
    
    def get_ENV(self,ENVFILENAME=None):
        hdf_env = h5py.File(ENVFILENAME)
        temperature = hdf_env['NS']['VERENV']['airTemperature'][:,12:37,:]
        da = xr.DataArray(temperature, dims=['along_track', 'cross_track','range'],
               coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                       'lats': (['along_track','cross_track'],self.xrds.lons),
                       'time': (['along_track','cross_track'],self.xrds.time),
                       'alt':(['along_track', 'cross_track','range'],self.xrds.alt)})

        da.where(da > 0)
        da.attrs['units'] =  'K'
        da.attrs['standard_name'] = 'GPM-DPR ENV data'
        self.xrds["T"] = da
    def get_merra(self,interp1=True,interp2=False,getsurf=True):
        """
        This method matches up the *closest* MERRA-2 profiles. 
        To do so it uses the xarray.sel command. 
        
        Please note this is not generalized. The files structure of my MERRA-2 files is a bit particular. 
        In theory you could point this into your own directory where those files are. Or even use a different
        reanalysis (e.g., ERA)
        
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

        ds_url = '/data/gpm/a/randyjc2/MERRA/PROFILE/'+ str(year) + '/' + 'MERRA2_400.inst6_3d_ana_Np.'+ str(year) + month + day+ '.nc4'

        ###load file
        merra = xr.open_dataset(ds_url,chunks={'lat': 361, 'lon': 576})
        ###

        #select the closest profile to the lat, lon, time
        sounding = merra.sel(lon=self.xrds.lons,lat=self.xrds.lats,time=self.xrds.time,method='nearest')
        sounding.load() 
        self.sounding = sounding
        
        if interp1:
            self.interp_MERRA(keyname='T')
            self.interp_MERRA(keyname='U')
            self.interp_MERRA(keyname='V')
            self.interp_MERRA(keyname='QV')
            self.interp_flag = 1
        elif interp2:
            self.interp_MERRA_V2(keyname='T')
            self.interp_MERRA_V2(keyname='U')
            self.interp_MERRA_V2(keyname='V')
            self.interp_MERRA_V2(keyname='QV')
            self.interp_flag = 1
            
        if getsurf:
            ds_url ='/data/keeling/a/randyjc2/gpm/MERRA/SURF/'+ str(year) + '/' + 'MERRA2_400.tavg1_2d_slv_Nx.'+str(year) + month + day +'.nc4'

            ###load file
            merra = xr.open_dataset(ds_url)
            ###

            #select the closest profile to the lat, lon, time
            gpmcoords = merra.sel(lon=self.xrds.lons,lat=self.xrds.lats,time=self.xrds.time,method='nearest')
            da = xr.DataArray(gpmcoords.T2M.values, dims=['along_track', 'cross_track'],
               coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                       'lats': (['along_track','cross_track'],self.xrds.lons),
                       'time': (['along_track','cross_track'],self.xrds.time)})
            da.attrs['units'] = gpmcoords.T2M.units
            da.attrs['standard_name'] = gpmcoords.T2M.standard_name
            self.xrds['T2M'] = da
            
            da = xr.DataArray(gpmcoords.T2MWET.values, dims=['along_track', 'cross_track'],
               coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                       'lats': (['along_track','cross_track'],self.xrds.lons),
                       'time': (['along_track','cross_track'],self.xrds.time)})
            da.attrs['units'] = gpmcoords.T2MWET.units
            da.attrs['standard_name'] = gpmcoords.T2MWET.standard_name
            self.xrds['T2MWET'] = da
            
            da = xr.DataArray(gpmcoords.T2MDEW.values, dims=['along_track', 'cross_track'],
               coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                       'lats': (['along_track','cross_track'],self.xrds.lons),
                       'time': (['along_track','cross_track'],self.xrds.time)})
            da.attrs['units'] = gpmcoords.T2MDEW.units
            da.attrs['standard_name'] = gpmcoords.T2MDEW.standard_name
            self.xrds['T2MDEW'] = da
            
            if self.snow:
                self.xrds = self.xrds.where(self.xrds.T2MWET-273.15 <= 0)
                #to reduce size of data, drop empty cross-track sections 
                self.xrds = self.xrds.dropna(dim='along_track',how='all')
                
            
    def interp_MERRA(self,keyname=None):
        """ 
        This interpolates the MERRA data from the self.get_merra method, to the same veritcal levels as the GPM-DPR
        
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
        return da
        
    def interp_MERRA_V2(self,keyname=None):
        """This is an effcient way of doing linear interpolation of the MERRA soundings """
        x = self.sounding['H'].values
        y = self.sounding[keyname].values
        z = self.xrds.alt.values*1000 #convert to m 
        interped = np.zeros(self.xrds.alt.values.shape)
        for i in np.arange(0,25):
            interped[:,i,:] = interp_2(x[:,i,:],y[:,i,:],z[0,i,:])

        da = xr.DataArray(interped, dims=['along_track', 'cross_track','range'],
               coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                       'lats': (['along_track','cross_track'],self.xrds.lons),
                       'time': (['along_track','cross_track'],self.xrds.time),
                       'alt':(['along_track', 'cross_track','range'],self.xrds.alt)})

        da.attrs['units'] = self.sounding[keyname].units
        da.attrs['standard_name'] = 'Interpolated ' + self.sounding[keyname].standard_name + ' to GPM height coord'
        self.xrds[keyname] = da
        
        return da
    
    def extract_nearsurf(self):
        """
        Since we are often concerned with whats happening at the surface, this will extract the variables just above
        the clutter. 
        """
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
            self.grab_variable(keyname='IWC')

#         if self.interp_flag == 1:
#             self.grab_variable(keyname='T')
#             self.grab_variable(keyname='U')
#             self.grab_variable(keyname='V')
#             self.grab_variable(keyname='QV')
        

    def grab_variable(self,keyname=None):
        """
        This goes along with the self.extract_nearsurf()
        """

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
                
    def get_physcial_distance(self,reference_point = None):
        """ 
        This method uses pyproj to calcualte distances between lats and lons. 
        reference_point is a list or array conisting of two entries, [Longitude,Latitude]
        
        Please note that this intentionally uses an older version of pyproj (< version 2.0, i used 1.9.5.1)
        This is because it preserves how the function is called. 
        """

        if reference_point is None and self.reference_point is None:
            print('Error, no reference point found...please enter one')
        else:
            #this is envoke the pyproj package. Please note this must be an old version** < 2.0 
            from pyproj import Proj
            p = Proj(proj='aeqd', ellps='WGS84', datum='WGS84', lat_0=reference_point[1], lon_0=reference_point[0])
            #double check to make sure this returns 0 meters
            x,y = p(reference_point[0],reference_point[1])
            if np.sqrt(x**2 + y**2) != 0:
                'something isnt right with the projection. investigate'
            else:
                ind = np.isnan(self.xrds.NSKu_nearSurf.values)
                x = np.zeros(self.xrds.lons.values.shape)
                y = np.zeros(self.xrds.lats.values.shape)
                x[~ind],y[~ind] = p(self.xrds.lons.values[~ind],self.xrds.lats.values[~ind])
                x[ind] = np.nan
                y[ind] = np.nan
                da = xr.DataArray(np.sqrt(x**2 + y**2)/1000, dims=['along_track', 'cross_track'],
                    coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                    'lats': (['along_track','cross_track'],self.xrds.lons),
                    'time': (['along_track','cross_track'],self.xrds.time)})
                da.attrs['units'] = 'km'
                da.attrs['standard_name'] = 'distance, way of the crow (i.e. direct), to the reference point'
                self.xrds['distance'] = da
                

def interp_2(x, y, z):
    """ This is from this discussion: https://stackoverflow.com/questions/14559687/scipy-fast-1-d-interpolation-without-any-loop"""
    rows, cols = x.shape
    row_idx = np.arange(rows).reshape((rows,) + (1,) * z.ndim)
    col_idx = np.argmax(x.reshape(x.shape + (1,) * z.ndim) > z, axis=1) - 1
    ret = y[row_idx, col_idx + 1] - y[row_idx, col_idx]
    ret /= x[row_idx, col_idx + 1] - x[row_idx, col_idx]
    ret *= z - x[row_idx, col_idx]
    ret += y[row_idx, col_idx]
    return ret
             
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
                
class APR():
    
    """
    Author: Randy J. Chase. This class is intended to help with APR-2/APR-3 files. 
    
    Currently supported campaigns: gcpex, olympex 
    
    Feel free to reach out to me on twitter (@dopplerchase) or email randyjc2@illinois.edu
    """
    
    def __init__(self):
        self.initialzied = True 
        self.T3d = False  
        
    def read(self,filename,campaign='gcpex'):
    
        """
        ===========

        This is for reading in apr3 hdf files from OLYMPEX and return them all in one dictionary

        ===========

        filename = filename of the apr3 file
        """
        
        if campaign=='gcpex':
            self.campaign = campaign
            from pyhdf.SD import SD, SDC
            apr = {}
            flag = 0
            ##Radar varibles in hdf file found by hdf.datasets
            radar_freq = 'zhh14' #Ku
            radar_freq2 = 'zhh35' #Ka
            radar_freq3 = 'zhh95' #W
            radar_freq4 = 'ldr14' #LDR
            vel_str = 'vel14' #Doppler
            ##

            hdf = SD(filename, SDC.READ)

            #What used to be here is a searching param looking for W-band, since we know there is no W in GCPEX, just autofill it. 
            alt = hdf.select('alt3D')
            lat = hdf.select('lat')
            lon = hdf.select('lon')
            roll = hdf.select('roll').get()
            time = hdf.select('scantime').get()
            surf = hdf.select('surface_index').get()
            isurf = hdf.select('isurf').get()
            plane = hdf.select('alt_nav').get()
            radar = hdf.select(radar_freq) #ku
            radar2 = hdf.select(radar_freq2) #ka
            radar4 = hdf.select(radar_freq4) #ldr
            vel = hdf.select(vel_str)
            lon3d = hdf.select('lon3D')
            lat3d = hdf.select('lat3D')
            alt3d = hdf.select('alt3D')
            lat3d_scale = hdf.select('lat3D_scale').get()[0][0]
            lon3d_scale = hdf.select('lon3D_scale').get()[0][0]
            alt3d_scale = hdf.select('alt3D_scale').get()[0][0]
            lat3d_offset = hdf.select('lat3D_offset').get()[0][0]
            lon3d_offset = hdf.select('lon3D_offset').get()[0][0]
            alt3d_offset = hdf.select('alt3D_offset').get()[0][0]

            alt = alt.get()
            ngates = alt.shape[0]
            lat = lat.get()
            lon = lon.get()

            lat3d = lat3d.get()
            lat3d = (lat3d/lat3d_scale) + lat3d_offset
            lon3d = lon3d.get()
            lon3d = (lon3d/lon3d_scale) + lon3d_offset
            alt3d = alt3d.get()
            alt3d = (alt3d/alt3d_scale) + alt3d_offset

            radar_n = radar.get()
            radar_n = radar_n/100.
            radar_n2 = radar2.get()
            radar_n2 = radar_n2/100.

            radar_n4 = radar4.get()
            radar_n4 = radar_n4/100.
            vel_n = vel.get()
            vel_n = vel_n/100.

            #Quality control (use masked values
            radar_n = np.ma.masked_where(radar_n<=-99,radar_n)
            radar_n2 = np.ma.masked_where(radar_n2<=-99,radar_n2)
            radar_n4 = np.ma.masked_where(radar_n4<=-99,radar_n4)          
            vel = np.ma.masked_where(vel_n <= -99,vel_n)

            #Get rid of nans, the new HDF has builtin
            radar_n = np.ma.masked_where(np.isnan(radar_n),radar_n)
            radar_n2 = np.ma.masked_where(np.isnan(radar_n2),radar_n2)
            radar_n4 = np.ma.masked_where(np.isnan(radar_n4),radar_n4)
            vel = np.ma.masked_where(np.isnan(vel),vel)

            radar_n3 = np.ma.zeros(radar_n.shape)
            radar_n3 = np.ma.masked_where(radar_n3 == 0,radar_n3)
        
        if campaign == 'olympex':
            
            self.campaign = campaign
            ##Radar varibles in hdf file found by hdf.datasets
            radar_freq = 'zhh14' #Ku
            radar_freq2 = 'zhh35' #Ka
            radar_freq3 = 'z95s' #W
            radar_freq4 = 'ldr14' #LDR
            vel_str = 'vel14' #Doppler
            ##


            import h5py
            hdf = h5py.File(filename,"r")

            listofkeys = hdf['lores'].keys()
            alt = hdf['lores']['alt3D'][:]
            lat = hdf['lores']['lat'][:]
            lon = hdf['lores']['lon'][:]
            time = hdf['lores']['scantime'][:]
            surf = hdf['lores']['surface_index'][:]
            isurf =  hdf['lores']['isurf'][:]
            plane =  hdf['lores']['alt_nav'][:]
            radar = hdf['lores'][radar_freq][:]
            radar2 = hdf['lores'][radar_freq2][:]
            radar4 = hdf['lores'][radar_freq4][:]
            vel = hdf['lores']['vel14c'][:]
            lon3d = hdf['lores']['lon3D'][:]
            lat3d = hdf['lores']['lat3D'][:]
            alt3d = hdf['lores']['alt3D'][:]
            roll = hdf['lores']['roll'][:]

            #see if there is W band
            if 'z95s' in listofkeys:
                if 'z95n' in listofkeys:
                    radar_nadir = hdf['lores']['z95n'][:]
                    radar_scanning = hdf['lores']['z95s'][:]
                    radar3 = radar_scanning
                    w_flag = 1
                    ##uncomment if you want high sensativty as nadir scan (WARNING, CALIBRATION)
                    #radar3[:,12,:] = radar_nadir[:,12,:]
                else:
                    radar3 = hdf['lores']['z95s'][:]
                    print('No vv, using hh')
            else:
                radar3 = np.ma.zeros(radar.shape)
                radar3 = np.ma.masked_where(radar3==0,radar3)
                w_flag = 1
                print('No W band')
            #Quality control (masked where invalid)
            radar_n = np.ma.masked_where(radar <= -99,radar)
            radar_n2 = np.ma.masked_where(radar2 <= -99,radar2)
            radar_n3 = np.ma.masked_where(radar3 <= -99,radar3)
            radar_n4 = np.ma.masked_where(radar4 <= -99,radar4)
            vel = np.ma.masked_where(vel <= -99,vel)

            #Get rid of nans, the new HDF has builtin
            radar_n = np.ma.masked_where(np.isnan(radar_n),radar_n)
            radar_n2 = np.ma.masked_where(np.isnan(radar_n2),radar_n2)
            radar_n3 = np.ma.masked_where(np.isnan(radar_n3),radar_n3)
            radar_n4 = np.ma.masked_where(np.isnan(radar_n4),radar_n4)
            vel = np.ma.masked_where(np.isnan(vel),vel)

        if campaign == 'camp2ex':
            ##Radar varibles in hdf file found by hdf.datasets
            radar_freq = 'zhh14' #Ku
            radar_freq2 = 'zhh35' #Ka
            radar_freq3 = 'zhh95' #W
            radar_freq4 = 'ldrhh14' #LDR
            vel_str = 'vel14' #Doppler
            ##
            
            import h5py
            hdf = h5py.File(filename,"r")

            listofkeys = hdf['lores'].keys()
            alt = hdf['lores']['alt3D'][:]
            lat = hdf['lores']['lat'][:]
            lon = hdf['lores']['lon'][:]
            time = hdf['lores']['scantime'][:]
            surf = hdf['lores']['surface_index'][:]
            isurf =  hdf['lores']['isurf'][:]
            plane =  hdf['lores']['alt_nav'][:]
            radar = hdf['lores'][radar_freq][:]
            radar2 = hdf['lores'][radar_freq2][:]
            radar4 = hdf['lores'][radar_freq4][:]
            vel = hdf['lores'][vel_str][:]
            lon3d = hdf['lores']['lon3D'][:]
            lat3d = hdf['lores']['lat3D'][:]
            alt3d = hdf['lores']['alt3D'][:]
            roll = hdf['lores']['roll'][:]

            #see if there is W band
            if 'z95s' in listofkeys:
                if 'z95n' in listofkeys:
                    radar_nadir = hdf['lores']['z95n'][:]
                    radar_scanning = hdf['lores']['z95s'][:]
                    radar3 = radar_scanning
                    w_flag = 1
                    ##uncomment if you want high sensativty as nadir scan (WARNING, CALIBRATION)
                    #radar3[:,12,:] = radar_nadir[:,12,:]
                else:
                    radar3 = hdf['lores']['z95s'][:]
                    print('No vv, using hh')
            else:
                radar3 = np.ma.zeros(radar.shape)
                radar3 = np.ma.masked_where(radar3==0,radar3)
                w_flag = 1
                print('No W band')
      

            #Quality control (masked where invalid)
            radar_n = np.ma.masked_where(radar <= -99,radar)
            radar_n2 = np.ma.masked_where(radar2 <= -99,radar2)
            radar_n3 = np.ma.masked_where(radar3 <= -99,radar3)
            radar_n4 = np.ma.masked_where(radar4 <= -99,radar4)
            vel = np.ma.masked_where(vel <= -99,vel)

            #Get rid of nans, the new HDF has builtin
            radar_n = np.ma.masked_where(np.isnan(radar_n),radar_n)
            radar_n2 = np.ma.masked_where(np.isnan(radar_n2),radar_n2)
            radar_n3 = np.ma.masked_where(np.isnan(radar_n3),radar_n3)
            radar_n4 = np.ma.masked_where(np.isnan(radar_n4),radar_n4)
            vel = np.ma.masked_where(np.isnan(vel),vel)
            
            
        ##convert time to datetimes
        time_dates = np.empty(time.shape,dtype=object)
        for i in np.arange(0,time.shape[0]):
            for j in np.arange(0,time.shape[1]):
                tmp = datetime.datetime.utcfromtimestamp(time[i,j])
                time_dates[i,j] = tmp

        #Create a time at each gate (assuming it is the same down each ray, there is a better way to do this)      
        time_gate = np.empty(lat3d.shape,dtype=object)
        for k in np.arange(0,lat3d.shape[0]):
            for i in np.arange(0,time_dates.shape[0]):
                for j in np.arange(0,time_dates.shape[1]):
                    time_gate[k,i,j] = time_dates[i,j]        
        
        time3d = np.copy(time_gate)

        da = xr.DataArray(radar_n,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                       'along_track':np.arange(radar_n.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],lon3d),
                                  'lat3d':(['range','cross_track','along_track'],lat3d),
                                  'time3d': (['range','cross_track','along_track'],time3d),
                                  'alt3d':(['range','cross_track','along_track'],alt3d)})      
        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Reflectivity'
        
        #make xr dataset
        self.xrds = da.to_dataset(name = 'Ku')
        #
        
        da = xr.DataArray(radar_n2,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                       'along_track':np.arange(radar_n.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],lon3d),
                                  'lat3d':(['range','cross_track','along_track'],lat3d),
                                  'time3d': (['range','cross_track','along_track'],time3d),
                                  'alt3d':(['range','cross_track','along_track'],alt3d)})  
        
        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ka-band Reflectivity'
        
        #add to  xr dataset
        self.xrds['Ka'] = da
        
        da = xr.DataArray(vel,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                       'along_track':np.arange(radar_n.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],lon3d),
                                  'lat3d':(['range','cross_track','along_track'],lat3d),
                                  'time3d': (['range','cross_track','along_track'],time3d),
                                  'alt3d':(['range','cross_track','along_track'],alt3d)})  
        
        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Doppler Vel'
        
        #add to  xr dataset
        self.xrds['DopKu'] = da
        #
        
        da = xr.DataArray(radar_n3,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                       'along_track':np.arange(radar_n.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],lon3d),
                                  'lat3d':(['range','cross_track','along_track'],lat3d),
                                  'time3d': (['range','cross_track','along_track'],time3d),
                                  'alt3d':(['range','cross_track','along_track'],alt3d)})  
        
        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'W-band Reflectivity'
        
        #add to  xr dataset
        self.xrds['W'] = da
        #
        
        da = xr.DataArray(radar_n4,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                       'along_track':np.arange(radar_n.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],lon3d),
                                  'lat3d':(['range','cross_track','along_track'],lat3d),
                                  'time3d': (['range','cross_track','along_track'],time3d),
                                  'alt3d':(['range','cross_track','along_track'],alt3d)})  
        
        da.fillna(value=-9999)
        da.attrs['units'] = 'dB'
        da.attrs['standard_name'] = 'LDR at Ku-band '
        
        #add to  xr dataset
        self.xrds['LDR'] = da
        #
        
        da = xr.DataArray(roll,dims={'cross_track':np.arange(0,24),'along_track':np.arange(radar_n.shape[2])})
        da.attrs['units'] = 'degrees'
        da.attrs['standard_name'] = 'Left/Right Plane Roll'
        
        #add to  xr dataset
        self.xrds['Roll'] = da
        #
        
        
        da = xr.DataArray(isurf,dims={'cross_track':np.arange(0,24),'along_track':np.arange(radar_n.shape[2])})
        da.attrs['units'] = 'none'
        da.attrs['standard_name'] = 'index of best guess for location of surfce'
        
        #add to  xr dataset
        self.xrds['surf'] = da
        #
            
    def determine_ground_lon_lat(self,near_surf_Z=True):
        
        mean_alt = self.xrds.alt3d.mean(axis=(1,2))
        ind_sealevel = find_nearest(mean_alt,0)
        g_lon = self.xrds.lon3d[ind_sealevel,:,:]
        g_lat = self.xrds.lat3d[ind_sealevel,:,:]
        
        da = xr.DataArray(g_lon,dims={'cross_track':g_lon.shape[0],'along_track':np.arange(g_lon.shape[1])})
        da.attrs['units'] = 'degress longitude'
        da.attrs['standard_name'] = 'Surface Longitude'
        self.xrds['g_lon'] = da

        da = xr.DataArray(g_lat,dims={'cross_track':g_lat.shape[0],'along_track':np.arange(g_lat.shape[1])})
        da.attrs['units'] = 'degress latitude'
        da.attrs['standard_name'] = 'Surface latitude'
        self.xrds['g_lat'] = da
        
        
        if near_surf_Z:
            ind_nearsurf = find_nearest(mean_alt,1100)
            g_Ku = self.xrds.Ku[ind_nearsurf,:,:]
            da = xr.DataArray(g_Ku,dims={'cross_track':g_lon.shape[0],'along_track':np.arange(g_lon.shape[1])})
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'Near_surf Z'
            self.xrds['Ku_nearsurf'] = da

            g_Ka = self.xrds.Ka[ind_nearsurf,:,:]
            da = xr.DataArray(g_Ka,dims={'cross_track':g_lon.shape[0],'along_track':np.arange(g_lon.shape[1])})
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'Near_surf Z'
            self.xrds['Ka_nearsurf'] = da

        
    def setboxcoords(self,nanout=True):
        """
        This method sets all points outside the box to nan. 
        """
        
        if len(self.corners) > 0:
            self.ll_lon = self.corners[0]
            self.ur_lon = self.corners[1]
            self.ll_lat = self.corners[2]
            self.ur_lat = self.corners[3]
            if nanout:
                self.xrds = self.xrds.where((self.xrds.g_lon >= self.ll_lon) & (self.xrds.g_lon <= self.ur_lon) & (self.xrds.g_lat >= self.ll_lat)  & (self.xrds.g_lat <= self.ur_lat),drop=False)
        else:
            print('ERROR, not boxcoods set...did you mean to do this?')
            
    def cit_temp_to_apr(self,fl,time_inds,insidebox=True):
    
        self.T3d = True 
        cit_lon = fl['longitude']['data'][time_inds]
        cit_lat = fl['latitude']['data'][time_inds]
        cit_alt = fl['altitude']['data'][time_inds]
        cit_twc =  fl['twc']['data'][time_inds]
        cit_T = fl['temperature']['data'][time_inds]
        
        if insidebox:
            ind_inbox = np.where((cit_lon >= self.ll_lon) & (cit_lon <= self.ur_lon) & (cit_lat >= self.ll_lat)  & (cit_lat <= self.ur_lat))
        else:
            ind_inbox = np.arange(0,len(fl['temperature']['data'][time_inds]))
            
        bins = np.arange(0,6000,500)
        binind = np.digitize(cit_alt[ind_inbox],bins=bins)
        df = pd.DataFrame({'Temperature':cit_T[ind_inbox],'Alt':cit_alt[ind_inbox],'binind':binind})
        df = df.groupby('binind').mean()
        f_T = scipy.interpolate.interp1d(df.Alt.values,df.Temperature.values,fill_value='extrapolate',kind='linear')
        T3d = f_T(self.xrds.alt3d.values)

        da = xr.DataArray(T3d,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})
        da.fillna(value=-9999)
        da.attrs['units'] = 'degC'
        da.attrs['standard_name'] = 'Temperature, inferred from Citation Spiral'

        #add to  xr dataset
        self.xrds['T3d'] = da
        #
       
    def correct_gaseous(self,filepathtogasploutput=None):
        """ This is a method to correct for 02 and H20 attenuation at Ku and Ka band, it requires you to run the gaspl package in matlab. If you wish to learn about this, please email me randyjc2 at illinois.edu """
        
        if filepathtogasploutput is None:
            print('Please supply filepath to gaspl output')
        
            return
        
        import scipy.io 
        import scipy.interpolate
        d = scipy.io.loadmat(filepathtogasploutput)
        #create interp funcs so we can plug in the apr gate structure
        ka_func = scipy.interpolate.interp1d(d['alt'].ravel(),d['L'].ravel(),kind='cubic',bounds_error=False) 
        ku_func = scipy.interpolate.interp1d(d['alt'].ravel(),d['L2'].ravel(),kind='cubic',bounds_error=False)
        
        k_ku = ku_func(self.xrds.alt3d.values)
        k_ka = ka_func(self.xrds.alt3d.values)

        k_ku = k_ku*0.03 #conver to db/gate
        k_ka = k_ka*0.03 #conver to db/gate

        k_ku[np.isnan(k_ku)] = 0
        k_ka[np.isnan(k_ka)] = 0

        k_ku = 2*np.cumsum(k_ku,axis=(0))
        k_ka = 2*np.cumsum(k_ka,axis=(0))
        
        ku_new = self.xrds.Ku.values + k_ku 
        da = xr.DataArray(ku_new,
                  dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                        'along_track':np.arange(self.xrds.Ku.shape[2])},
                  coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                          'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                          'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                          'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Reflectivity'

        self.xrds['Ku'] = da
        
        ka_new = self.xrds.Ka.values + k_ka
        da = xr.DataArray(ka_new,
                  dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                        'along_track':np.arange(self.xrds.Ku.shape[2])},
                  coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                          'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                          'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                          'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ka-band Reflectivity'

        self.xrds['Ka'] = da
        
        return
        
    def run_retrieval(self,old=True):
        
        from pickle import load
        if old:
            scaler_X = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/scaler_X.pkl', 'rb'))
            scaler_y = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/scaler_y.pkl', 'rb'))
        else:
            scaler_X = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/scaler_X_V2.pkl', 'rb'))
            scaler_y = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/scaler_y_V2.pkl', 'rb'))
            
        #now we have to reshape things to make sure they are in the right shape for the NN model [n_samples,n_features]
        Ku = self.xrds.Ku.values
        shape_step1 = Ku.shape
        Ku = Ku.reshape([Ku.shape[0],Ku.shape[1]*Ku.shape[2]])
        shape_step2 = Ku.shape
        Ku = Ku.reshape([Ku.shape[0]*Ku.shape[1]])
        Ka = self.xrds.Ka.values
        Ka = Ka.reshape([Ka.shape[0],Ka.shape[1]*Ka.shape[2]])
        Ka = Ka.reshape([Ka.shape[0]*Ka.shape[1]])

        T = self.xrds.T3d.values
        T = T.reshape([T.shape[0],T.shape[1]*T.shape[2]])
        T = T.reshape([T.shape[0]*T.shape[1]])

        ind_masked = np.isnan(Ku) 
        Ku_nomask = np.zeros(Ku.shape)
        Ka_nomask = np.zeros(Ka.shape)
        T_nomask = np.zeros(T.shape)

        Ku_nomask[~ind_masked] = Ku[~ind_masked]
        Ka_nomask[~ind_masked] = Ka[~ind_masked]
        T_nomask[~ind_masked]  = T[~ind_masked]
        ind = np.where(Ku_nomask!=0)[0]

        #scale the input vectors by the mean that it was trained with
        X = np.zeros([Ku_nomask.shape[0],3])
        X[:,0] = (Ku_nomask - scaler_X.mean_[0])/scaler_X.scale_[0] #ku 
        X[:,1] = ((Ku_nomask - Ka_nomask)- scaler_X.mean_[1])/scaler_X.scale_[1] #dfr
        X[:,2] = (T_nomask - scaler_X.mean_[2])/scaler_X.scale_[2] #T
        #

        import tensorflow as tf
        from tensorflow.python.keras import losses
        if old:
            model=tf.keras.models.load_model('/data/gpm/a/randyjc2/DRpy/drpy/models/NN_4by8.h5',
                                     custom_objects=None,compile=True)
        else:
            model=tf.keras.models.load_model('/data/gpm/a/randyjc2/DRpy/drpy/models/NN_6by8.h5',
                                     custom_objects=None,compile=True)
        if self.T3d:
            yhat = model.predict(X[ind,0:3],batch_size=len(X[ind,0]))
            yhat = scaler_y.inverse_transform(yhat)
            yhat[:,1] = 10**yhat[:,1]
            yhat[:,2] = 10**yhat[:,2]
        else:
            print('ERROR, No temperautre data')
        
        ind = np.where(Ku_nomask!=0)[0]
        Nw = np.zeros(Ku_nomask.shape)
        Nw[ind] = np.squeeze(yhat[:,0])
        Nw = Nw.reshape([shape_step2[0],shape_step2[1]])
        Nw = Nw.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
        Nw = np.ma.masked_where(Nw==0.0,Nw)

        da = xr.DataArray(Nw,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})
        da.fillna(value=-9999)
        da.attrs['units'] = 'log(m^-4)'
        da.attrs['standard_name'] = 'Retireved Liquid Eq. Nw'

        self.xrds['Nw'] = da

        Dm = np.zeros(Ku_nomask.shape)
        Dm[ind] = np.squeeze(yhat[:,1])
        Dm = Dm.reshape([shape_step2[0],shape_step2[1]])
        Dm = Dm.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
        Dm = np.ma.masked_where(Dm==0.0,Dm)

        da = xr.DataArray(Dm,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'mm'
        da.attrs['standard_name'] = 'Retireved Liquid Eq. Dm'

        self.xrds['Dm'] = da
        
        Dm_frozen = np.zeros(Ku_nomask.shape)
        Dm_frozen[ind] = np.squeeze(yhat[:,2])
        Dm_frozen = Dm_frozen.reshape([shape_step2[0],shape_step2[1]])
        Dm_frozen = Dm_frozen.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
        Dm_frozen = np.ma.masked_where(Dm_frozen==0.0,Dm_frozen)

        da = xr.DataArray(Dm_frozen,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'mm'
        da.attrs['standard_name'] = 'Retireved Frozen Dm'

        self.xrds['Dm_frozen'] = da
        
        Nw = 10**Nw #undo log, should be in m^-4
        Dm = Dm/1000. # convert to m ^4
        IWC = (Nw*(Dm)**4*1000*np.pi)/4**(4) # the 1000 is density of water (kg/m^3)
        IWC = IWC*1000 #convert to g/m^3 
        
        da = xr.DataArray(IWC,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'g m^-3'
        da.attrs['standard_name'] = 'Retireved IWC, calculated from Nw and Dm'

        self.xrds['IWC'] = da
        
    def get_legacy(self):
        
        from pickle import load
        #load GPM-DPR v06 average relations 
        reg_kupr = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/LinearRegression_KUPR.pkl', 'rb'))
        reg_kapr = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/LinearRegression_KAPR.pkl', 'rb'))
        #load newly trained exp relations
        reg_ku = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/LinearRegression_KU.pkl', 'rb'))
        reg_ka = load(open('/data/gpm/a/randyjc2/DRpy/drpy/models/LinearRegression_KA.pkl', 'rb'))

        #now we have to reshape things to make sure they are in the right shape for the NN model [n_samples,n_features]
        Ku = 10**(self.xrds.Ku.values/10)
        shape_step1 = Ku.shape
        Ku = Ku.reshape([Ku.shape[0],Ku.shape[1]*Ku.shape[2]])
        shape_step2 = Ku.shape
        Ku = Ku.reshape([Ku.shape[0]*Ku.shape[1]])
        Ka = 10**(self.xrds.Ka.values/10)
        Ka = Ka.reshape([Ka.shape[0],Ka.shape[1]*Ka.shape[2]])
        Ka = Ka.reshape([Ka.shape[0]*Ka.shape[1]])

        ind_masked = np.isnan(Ku) 
        Ku_nomask = np.zeros(Ku.shape)
        Ka_nomask = np.zeros(Ka.shape)

        Ku_nomask[~ind_masked] = Ku[~ind_masked]
        Ka_nomask[~ind_masked] = Ka[~ind_masked]
        ind = np.where(Ku_nomask!=0)[0]

        #Do Ku-band relations 
        X = np.zeros([Ku_nomask.shape[0],1])
        X[:,0] = Ku_nomask

        #kuexp 
        yhat = 10**reg_ku.predict(np.log10(X[ind,:]))

        IWC = np.zeros(Ku_nomask.shape)
        IWC[ind] = np.squeeze(yhat)
        IWC = IWC.reshape([shape_step2[0],shape_step2[1]])
        IWC = IWC.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
        IWC = np.ma.masked_where(IWC==0.0,IWC)

        da = xr.DataArray(IWC,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'g m^-3'
        da.attrs['standard_name'] = 'Retireved IWC from new exp (ku-band)'

        self.xrds['IWC_kuexp'] = da

        #kuprexp 
        yhat = 10**reg_kupr.predict(np.log10(X[ind,:]))

        IWC = np.zeros(Ku_nomask.shape)
        IWC[ind] = np.squeeze(yhat)
        IWC = IWC.reshape([shape_step2[0],shape_step2[1]])
        IWC = IWC.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
        IWC = np.ma.masked_where(IWC==0.0,IWC)

        da = xr.DataArray(IWC,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'g m^-3'
        da.attrs['standard_name'] = 'Retireved IWC from kupr exp'

        self.xrds['IWC_kupr'] = da

        #Do Ka-band relations 
        X = np.zeros([Ku_nomask.shape[0],1])
        X[:,0] = Ka_nomask

        #kuexp 
        yhat = 10**reg_ka.predict(np.log10(X[ind,:]))

        IWC = np.zeros(Ku_nomask.shape)
        IWC[ind] = np.squeeze(yhat)
        IWC = IWC.reshape([shape_step2[0],shape_step2[1]])
        IWC = IWC.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
        IWC = np.ma.masked_where(IWC==0.0,IWC)

        da = xr.DataArray(IWC,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'g m^-3'
        da.attrs['standard_name'] = 'Retireved IWC from new exp (ka-band)'

        self.xrds['IWC_kaexp'] = da

        #kuprexp 
        yhat = 10**reg_kapr.predict(np.log10(X[ind,:]))

        IWC = np.zeros(Ku_nomask.shape)
        IWC[ind] = np.squeeze(yhat)
        IWC = IWC.reshape([shape_step2[0],shape_step2[1]])
        IWC = IWC.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
        IWC = np.ma.masked_where(IWC==0.0,IWC)

        da = xr.DataArray(IWC,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'g m^-3'
        da.attrs['standard_name'] = 'Retireved IWC from kapr exp'

        self.xrds['IWC_kapr'] = da
    def mask_lobe(self):
        """ MASK SIDELOBES. This is a problem in GCPEX. NEED TO MAKE SURE THEY ARE GONE """ 
        ku_new = np.ma.masked_where(np.ma.abs(self.xrds.DopKu.values) >= 10, self.xrds.Ku.values)
        da = xr.DataArray(ku_new,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Reflectivity'

        self.xrds['Ku'] = da
        
    def cloudtopmask(self,sigma=1,mindBZ = 10,mask_others=True):
        import scipy.ndimage
        ku_temp = self.xrds.Ku.values 
        ku_temp[np.isnan(ku_temp)] = -99.99
        #Maskes Ku cloudtop noise
        ku_new = np.zeros(ku_temp.shape)
        ku_new = np.ma.masked_where(ku_new == 0.0,ku_new)
        for i in np.arange(0,24):
            temp = np.copy(ku_temp[:,i,:])
            a = scipy.ndimage.filters.gaussian_filter(temp,sigma)
            a = np.ma.masked_where(a<mindBZ,a)
            temp = np.ma.masked_where(a.mask,temp)
            ku_new[:,i,:] = np.ma.copy(temp)

        # np.ma.set_fill_value(ku_temp,-9999)
        ku_new = np.ma.masked_where(ku_new < mindBZ,ku_new)
        #new data has some weird data near plane. Delete manually
        ku_new = np.ma.masked_where(self.xrds.alt3d.values > 10000, ku_new)
        
        da = xr.DataArray(ku_new,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Reflectivity'

        self.xrds['Ku'] = da
        
        if mask_others:
            ka_new = np.ma.masked_where(ku_new.mask,self.xrds.Ka.values)
            da = xr.DataArray(ka_new,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

            da.fillna(value=-9999)
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'Ka-band Reflectivity'

            self.xrds['Ka'] = da
            
            w_new = np.ma.masked_where(ku_new.mask,self.xrds.W.values)
            da = xr.DataArray(w_new,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

            da.fillna(value=-9999)
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'W-band Reflectivity'

            self.xrds['W'] = da

    def rollfix(self):
            roll = self.xrds.Roll.values
            roll3d = roll[np.newaxis,:,:]
            roll3d = np.tile(roll3d,(self.xrds.Ku.values.shape[0],1,1))
            ku_new = np.ma.masked_where(np.abs(roll3d) > 10, self.xrds.Ku.values)
            da = xr.DataArray(ku_new,
                  dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                        'along_track':np.arange(self.xrds.Ku.shape[2])},
                  coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                          'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                          'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                          'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})
            da.fillna(value=-9999)
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'Ku-band Reflectivity'
            ka_new = np.ma.masked_where(ku_new.mask,self.xrds.Ka.values)
            da = xr.DataArray(ka_new,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})
            da.fillna(value=-9999)
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'Ka-band Reflectivity'
            self.xrds['Ka'] = da
            w_new = np.ma.masked_where(ku_new.mask,self.xrds.W.values)
            da = xr.DataArray(w_new,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

            da.fillna(value=-9999)
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'W-band Reflectivity'
            self.xrds['W'] = da
            
    def convolve_ku(self,x_stddev=1,y_stddev=1,x_size=3,y_size=3):
        
        from astropy.convolution import convolve
        from astropy.convolution import Gaussian2DKernel
        
        kernel = Gaussian2DKernel(x_stddev=x_stddev,y_stddev=y_stddev,x_size=x_size,y_size=y_size)
        ku = self.xrds.Ku.values
        ku_new = np.zeros(ku.shape)
        for i in np.arange(0,23):
            ku_new[:,i,:] = convolve(ku[:,i,:].squeeze(),kernel)
        ku_new = np.ma.masked_where(self.xrds.alt3d.values > 10000, ku_new)
        da = xr.DataArray(ku_new,
                  dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                        'along_track':np.arange(self.xrds.Ku.shape[2])},
                  coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                          'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                          'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                          'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Reflectivity'
        self.xrds['Ku'] = da
        
    def remove_skin_paint(self,n_convolutions=2,x_stddev=1,y_stddev=1,x_size=3,y_size=3):
        """This is a method to remove aircraft skinpaints. LDR shouldn't be avail for any good snow echoes"""
        from astropy.convolution import convolve
        from astropy.convolution import Gaussian2DKernel

        kernel = Gaussian2DKernel(x_stddev=x_stddev,y_stddev=y_stddev,x_size=x_size,y_size=y_size)
        
        ldr = self.xrds.LDR.values
        for j in np.arange(0,n_convolutions):
            ldr_new = np.zeros(ldr.shape)
            for i in np.arange(0,23):
                ldr_new[:,i,:] = convolve(ldr[:,i,:].squeeze(),kernel)

            ldr = np.copy(ldr_new)

        ku = self.xrds.Ku.values
        ku_new = np.ma.masked_where(~np.isnan(ldr_new),ku)
        da = xr.DataArray(ku_new,
                  dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                        'along_track':np.arange(self.xrds.Ku.shape[2])},
                  coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                          'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                          'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                          'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Reflectivity'
        self.xrds['Ku'] = da
        
    def remove_skin_paint_mannual(self,xind=[0,0],yind=[0,0],zind=[0,0]):
        """ indices are determined from near nadir """
        import scipy.interpolate
        ku = np.copy(self.xrds.Ku.values[:,:,:])
        for z in np.arange(zind[0],zind[1]):
            for y in np.arange(0,23):
                f = scipy.interpolate.interp1d(np.arange(0,550),
                                               self.xrds.alt3d.values[:,12,z+self.ind_start],
                                               bounds_error=False) #find the alt of the index at nadir
                f2 = scipy.interpolate.interp1d(self.xrds.alt3d.values[:,y,z+self.ind_start],
                                               np.arange(0,550),
                                               bounds_error=False) #find the index of the alt off nadir 
                alt_rm = f(xind[0])
                ind_adjust = f2(alt_rm)
                ind_low = np.asarray(np.round(ind_adjust),dtype=int)

                alt_rm = f(xind[1])
                ind_adjust = f2(alt_rm)
                ind_high = np.asarray(np.round(ind_adjust),dtype=int)
                ku[ind_low:ind_high,y,z+self.ind_start] = np.nan  

        da = xr.DataArray(ku,
                  dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                        'along_track':np.arange(self.xrds.Ku.shape[2])},
                  coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                          'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                          'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                          'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Reflectivity'
        self.xrds['Ku'] = da   
        
    def mask_surf(self,mask_others=True):
        k = 12
        ku_orig = self.xrds.Ku.values[:,k,:]
        ku_out = np.ma.zeros(self.xrds.Ku.shape)
        ku_orig = np.ma.masked_where(np.isnan(ku_orig),ku_orig)
        ku = scipy.ndimage.gaussian_filter(ku_orig,1)
        for i in np.arange(0,ku.shape[1]):
            ku_prof = 10**(ku[:,i]/10.)
            alt_prof = self.xrds.alt3d.values[:,k,i]
            surf1 = self.xrds.surf.values[k,i]
            filt = scipy.ndimage.sobel(ku_prof)
            #from the surface go up until Der < some value
            filt2 = np.abs(filt)
#             from IPython.core.debugger import Tracer; Tracer()() 
            tol = filt2[int(surf1)]/filt2[int(surf1)]
            j = 0 
            while tol > 0.01:
                j = j - 1
                ii = int(surf1)+j
                tol = filt2[ii]/filt2[int(surf1)]
            val = 2500
            if alt_prof[ii] >= val:
                ku_prof = np.ma.masked_where(alt_prof <= val,ku_orig[:,i])
                ku_orig[:,i] = ku_prof
            else:
                ku_orig[ii:,i] = np.ma.masked
        
        ku_out[:,k,:] = ku_orig

        for k in np.arange(0,24):
            if k == 12:
                continue
            ku_out[:,k,:] = np.ma.masked_where(ku_out[:,12,:].mask,self.xrds.Ku.values[:,k,:])
        
        if self.campaign=='gcpex':
            ku_out = np.ma.masked_where(self.xrds.alt3d.values < 600,ku_out)
            
        da = xr.DataArray(ku_out,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Reflectivity'
        self.xrds['Ku'] = da
        
        if mask_others:
            ka_new = np.ma.masked_where(ku_out.mask,self.xrds.Ka.values)
            da = xr.DataArray(ka_new,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

            da.fillna(value=-9999)
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'Ka-band Reflectivity'

            self.xrds['Ka'] = da
            
            w_new = np.ma.masked_where(ku_out.mask,self.xrds.W.values)
            da = xr.DataArray(w_new,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

            da.fillna(value=-9999)
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'W-band Reflectivity'

            self.xrds['W'] = da
    
    def mask_BB(self,mask_others=True):
        
        ldr = np.ma.masked_where(np.isnan(self.xrds.Ku.values),self.xrds.LDR.values)
        bb = precip_echo_filt3D(ldr,thresh=7)
        ind1 = np.where(bb[12,:] == 1) #BB profiles based on LDR
        top_a = find_bb(self.xrds,ind1)
        bb_long = extend_bb(ind1,self.xrds.time3d.values[0,12,:],top_a)
        #store bb_long for constraining Citation 
        self.bb_long = bb_long 
        ku_new  = np.ma.masked_where(self.xrds.alt3d.values <= bb_long,self.xrds.Ku.values)
        da = xr.DataArray(ku_new,
                  dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                        'along_track':np.arange(self.xrds.Ku.shape[2])},
                  coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                          'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                          'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                          'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})
        
        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Reflectivity'
        self.xrds['Ku'] = da
        
        if mask_others:
            ka_new = np.ma.masked_where(ku_new.mask,self.xrds.Ka.values)
            da = xr.DataArray(ka_new,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

            da.fillna(value=-9999)
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'Ka-band Reflectivity'

            self.xrds['Ka'] = da
            
            w_new = np.ma.masked_where(ku_new.mask,self.xrds.W.values)
            da = xr.DataArray(w_new,
                          dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                                'along_track':np.arange(self.xrds.Ku.shape[2])},
                          coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                                  'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                                  'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                                  'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

            da.fillna(value=-9999)
            da.attrs['units'] = 'dBZ'
            da.attrs['standard_name'] = 'W-band Reflectivity'

            self.xrds['W'] = da
         

    def mask_each_other(self):
        """ Make sure the Ku and Ka band radar are masked in the same spots """
        
        ku = np.copy(self.xrds.Ku.values)
        ka = np.copy(self.xrds.Ka.values)

        ka_new = np.ma.masked_where(np.isnan(ku),ka)
        ku_new = np.ma.masked_where(np.isnan(ka_new),ku)

        da = xr.DataArray(ka_new,
                      dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                            'along_track':np.arange(self.xrds.Ku.shape[2])},
                      coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                              'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                              'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                              'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ka-band Reflectivity'

        self.xrds['Ka'] = da

        da = xr.DataArray(ku_new,
                      dims={'range':np.arange(0,550),'cross_track':np.arange(0,24),
                            'along_track':np.arange(self.xrds.Ku.shape[2])},
                      coords={'lon3d': (['range','cross_track','along_track'],self.xrds.lon3d),
                              'lat3d': (['range','cross_track','along_track'],self.xrds.lat3d),
                              'time3d': (['range','cross_track','along_track'],self.xrds.time3d),
                              'alt3d':(['range','cross_track','along_track'],self.xrds.alt3d)})

        da.fillna(value=-9999)
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'Ku-band Reflectivity'

        self.xrds['Ku'] = da
            
def find_bb(ds,ind1):
    ku = np.squeeze(ds.Ku.values[:,12,ind1])
    alt = np.squeeze(ds.alt3d.values[:,12,ind1])
    surf = ds.surf.values[12,ind1]
    ind = np.isnan(ku) 
    ku[ind] = 0
    ku = scipy.ndimage.gaussian_filter(ku,5)
    top_a = np.zeros(ku.shape[1])
    for i in np.arange(0,ku.shape[1]):
        a = i
        ku_prof = 10**(ku[:,a]/10.)
        alt_prof = alt[:,a]

        filt = scipy.ndimage.sobel(ku_prof)
        filtfilt = scipy.ndimage.sobel(filt)
        filtfiltfilt = scipy.ndimage.sobel(filtfilt)
        k_func1 = scipy.interpolate.interp1d(alt_prof,ku_prof,kind='cubic',bounds_error=False)
        k_func2 = scipy.interpolate.interp1d(alt_prof,filtfilt,kind='cubic',bounds_error=False) 
        alt_new = np.linspace(alt_prof.min(),alt_prof.max(),1000)
        ku_interp = k_func1(alt_new)
        der2 = k_func2(alt_new)
        bbmax = np.where(ku_interp == ku_interp.max())

        ind = np.where(alt_new >= alt_new[bbmax])
        ind2 = np.where(alt_new <= alt_new[bbmax]+500)
        ind3 = np.intersect1d(ind,ind2)

        der2_sub = der2[ind3]
        alt_sub = alt_new[ind3]

        k_func3 = scipy.interpolate.interp1d(alt_sub,der2_sub,kind='cubic',bounds_error=False)
        try:
            top  = scipy.optimize.bisect(k_func3,alt_sub[0],alt_sub[-1])
        except:
            top  = scipy.optimize.bisect(k_func3,alt_sub[0],alt_sub[-1]+100)

        top_a[a] = top

    #clean up top 
    ind = np.where(np.abs(top_a-np.percentile(top_a,50)) >= 300)
    top_a[ind] = np.percentile(top_a,50)

    return top_a

def extend_bb(ind,R,bb_center):

    tots = np.arange(0,len(R))
    ex = np.setdiff1d(tots,ind)
    both = np.intersect1d(ind,tots)

    bb_long = np.ma.zeros(len(R))
    bb_long[ind] = bb_center
    for i in np.arange(0,len(ex)):
        t = ex[i]
        index = find_nearest(ind[0],t)
        bb_long[t] = bb_long[ind[0][index]]
    
    return bb_long
        
def precip_echo_filt3D(ku,thresh=5):
    """
    
    This function provides a filter for APR3 data to determine if there exists a precip echo in the column. Only preforms it 
    on one scan at a time. Could easily go up to all scans. 
    
    MAKE SURE TO APPLY masking routines first... 
    
    """
    precip_yn = np.zeros([ku.shape[1],ku.shape[2]])
    for k in np.arange(0,ku.shape[1]):
        for j in np.arange(0,ku.shape[2]):
            flag1 = 0 
            c1 = -1
            i = -1
            start_flag = 0
            while flag1 == 0:
                i = i + 1

                if c1 >= thresh:
                    precip_yn[k,j] = 1
                    break
                if i == 550: 
                    precip_yn[k,j] = 0
                    break

                t = ku[i,k,j]

                if start_flag ==0:
                    if np.ma.is_masked(t):
                        continue
                    else:
                        start_flag = 1
                        c1 = c1 + 1
                else:
                    if np.ma.is_masked(t):
                        start_flag = 0
                        c1 = -1
                        continue
                    else:
                        c1 = c1 + 1
    return precip_yn
        
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx