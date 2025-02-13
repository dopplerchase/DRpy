from __future__ import absolute_import
import xarray as xr 
import numpy as np
import pandas as pd
import datetime
import os
#turn off warnings so i can use the progressbar
import warnings
warnings.filterwarnings('ignore')

class GPMDPR():

    """
    Author: Randy J. Chase. This class is intended to help with the efficient processing of GPM-DPR radar files. 
    Xarray now supports direct reading the NASA hdf5 files. Here is the major refactor from the previous version 
    of DRpy. 

    Feel free to reach out to me on twitter (@dopplerchase) or email randychase@ou.edu
    
    For your reference, please check out GPM-DPR's ATBD: https://pps.gsfc.nasa.gov/GPMprelimdocs.html 
    """

    def __init__(self,filename=[],bounding_box=None,outer_swath=False,auto_run=True,heavy=True): 
        """
        Initializes things

        params::
        filename: str, path to GPM-DPR file 
        boundingbox: list of floats, if you would like to cut the gpm to a lat lon box 
        send in a list of [lon_min,lon_mat,lat_min,lat_max]
        """
        self.filename = filename
        self.corners = bounding_box
        self.heavy=heavy
        
        if auto_run:
            #this reads the hdf5 file 
            self.read()
            #this gets a datetime obj for each scan 
            self.parse_dtime()
        
    def read(self):
        """
        This method unfolds all the groups into one combined xarray dataset for
        each radar (e.g., KuPR and KaPR). It will be lazily loaded to save on RAM.
        Note that this code was primarily developed for V7 DPR products. It will not 
        work for V6 data 

        """
        #######################################################################
        ################################ KuPR #################################
        #######################################################################

        prefix = '/FS/'
        geo = xr.open_dataset(self.filename,group=prefix,engine='netcdf4',decode_cf=False)
        pre = xr.open_dataset(self.filename,group=prefix+'PRE',engine='netcdf4',decode_cf=False)
        slv = xr.open_dataset(self.filename,group=prefix+'SLV',engine='netcdf4',decode_cf=False)
        tim = xr.open_dataset(self.filename,group=prefix+'ScanTime',engine='netcdf4',decode_cf=False)

        #rename dims to proper names
        bad_dims = list(geo.dims)
        geo = geo.rename_dims({bad_dims[0]:'nscan',
                            bad_dims[1]:'nrayNS'})
        bad_dims = list(pre.dims)
        pre = pre.rename_dims({bad_dims[0]:'nscan',
                            bad_dims[1]:'nrayNS',
                            bad_dims[2]:'nfreq',
                            bad_dims[3]:'nbin'})
        bad_dims = list(slv.dims)
        slv = slv.rename_dims({bad_dims[0]:'nscan',
                            bad_dims[1]:'nrayNS',
                            bad_dims[2]:'nbin',
                            bad_dims[3]:'nfreq',
                            bad_dims[4]:'nNUBF'})
        bad_dims = list(tim.dims)
        tim = tim.rename_dims({bad_dims[0]:'nscan'})
        
        #if you want to load the whole thing in, turn the heavy flag on 
        if self.heavy: 
            ver = xr.open_dataset(self.filename,group=prefix+'VER',engine='netcdf4',decode_cf=False)
            srt = xr.open_dataset(self.filename,group=prefix+'SRT',engine='netcdf4',decode_cf=False)
            csf = xr.open_dataset(self.filename,group=prefix+'CSF',engine='netcdf4',decode_cf=False)
            exp = xr.open_dataset(self.filename,group=prefix+'Experimental',engine='netcdf4',decode_cf=False)
            flg = xr.open_dataset(self.filename,group=prefix+'FLG',engine='netcdf4',decode_cf=False)
            trg = xr.open_dataset(self.filename,group=prefix+'TRG',engine='netcdf4',decode_cf=False)

            #rename dims to proper names
            bad_dims = list(ver.dims)
            ver = ver.rename_dims({bad_dims[0]:'nscan',
                                bad_dims[1]:'nrayNS',
                                bad_dims[2]:'nbin',
                                bad_dims[3]:'nfreq',
                                bad_dims[4]:'nNP'})
            bad_dims = list(srt.dims)
            srt = srt.rename_dims({bad_dims[0]:'nscan',
                                bad_dims[1]:'nrayNS',
                                bad_dims[2]:'method',
                                bad_dims[3]:'foreBack',
                                bad_dims[4]:'nearFar',
                                bad_dims[5]:'nsdew'})
            bad_dims = list(csf.dims)
            csf = csf.rename_dims({bad_dims[0]:'nscan',
                                bad_dims[1]:'nrayNS',
                                bad_dims[2]:'nfreqHI'})
            bad_dims = list(exp.dims)
            exp = exp.rename_dims({bad_dims[0]:'nscan',
                                bad_dims[1]:'nrayNS',
                                bad_dims[2]:'nbinSZP',
                                bad_dims[3]:'nfreq'}) 
            bad_dims = list(flg.dims)
            flg = flg.rename_dims({bad_dims[0]:'nscan',
                                bad_dims[1]:'nrayNS',
                                bad_dims[2]:'nbin',
                                bad_dims[3]:'nfreq'})       

            bad_dims = list(trg.dims)
            trg = trg.rename_dims({bad_dims[0]:'nscan',
                                bad_dims[1]:'nrayNS',
                                bad_dims[2]:'nslope',})     

            #MERGE into one ds 
            self.ds = xr.merge([geo,pre,slv,ver,srt,csf,exp,flg,tim,trg])
        else:
            #MERGE into one ds 
            self.ds = xr.merge([geo,pre,slv,tim])
        

        #close uneeded xr datasets 
        geo.close()
        pre.close()
        slv.close()
        tim.close()
        if self.heavy: 
          ver.close()
          srt.close()
          csf.close()
          exp.close()
          flg.close()
        
        #set lat,lon,height as the coords to allow for easy xr slicing
        self.ds = self.ds.set_coords(['Latitude','Longitude','height'])

    def setboxcoords(self):
        """
        This method sets all points outside the box to nan. 
        """
        if len(self.corners) > 0:
            self.ll_lon = self.corners[0]
            self.ur_lon = self.corners[1]
            self.ll_lat = self.corners[2]
            self.ur_lat = self.corners[3]
            self.ds = self.ds.where((self.ds.Longitude >= self.ll_lon) & (self.ds.Longitude <= self.ur_lon) & (self.ds.Latitude >= self.ll_lat)  & (self.ds.Latitude <= self.ur_lat),drop=False)
        else:
            print('ERROR, no boxcoods set...did you mean to do this?')
    def parse_dtime(self):
        year = self.ds.Year.values
        ind = np.where(year == -9999)[0]
        year = np.asarray(year,dtype=str)
        year = list(year)

        month = self.ds.Month.values
        month = np.asarray(month,dtype=str)
        month = np.char.rjust(month, 2, fillchar='0')
        month = list(month)

        day = self.ds.DayOfMonth.values
        day = np.asarray(day,dtype=str)
        day = np.char.rjust(day, 2, fillchar='0')
        day = list(day)

        hour = self.ds.Hour.values
        hour = np.asarray(hour,dtype=str)
        hour = np.char.rjust(hour, 2, fillchar='0')
        hour = list(hour)

        minute = self.ds.Minute.values
        minute = np.asarray(minute,dtype=str)
        minute = np.char.rjust(minute, 2, fillchar='0')
        minute = list(minute)

        second = self.ds.Second.values
        second = np.asarray(second,dtype=str)
        second = np.char.rjust(second, 2, fillchar='0')
        second = list(second)

        millisecond = self.ds.MilliSecond.values
        millisecond = np.asarray(millisecond,dtype=str)
        millisecond = np.char.rjust(millisecond, 3, fillchar='0')
        millisecond = list(millisecond)

        datestr  = [year[i] +"-"+ month[i]+ "-" + day[i] + \
                    ' ' + hour[i] + ':' + minute[i] + ':' + second[i] + '.' + millisecond[i] for i in range(len(year))]
        datestr = np.asarray(datestr,dtype=str)
        datestr[ind] = '1970-01-01 00:00:00'
        datestr = np.reshape(datestr,[len(datestr),1])
        datestr = np.tile(datestr,(1,49))

        self.ds['time'] = xr.DataArray(np.asarray(datestr,dtype=np.datetime64),dims=['nscan','nrayNS'])
        self.ds = self.ds.set_coords('time')
        #drop the variables we dont need 
        self.ds = self.ds.drop(['Year','Month','DayOfMonth','Hour','Minute','Second','MilliSecond','DayOfYear','SecondOfDay'])

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
                ind = np.isnan(self.ds.zFactorFinalNearSurface.values[:,:,0])
                x = np.zeros(self.ds.Longitude.values.shape)
                y = np.zeros(self.ds.Latitude.values.shape)
                x[~ind],y[~ind] = p(self.ds.Longitude.values[~ind],self.ds.Latitude.values[~ind])
                x[ind] = np.nan
                y[ind] = np.nan
                da = xr.DataArray(np.sqrt(x**2 + y**2)/1000, dims=['nscan', 'nrayNS'])
                da.attrs['units'] = 'km'
                da.attrs['standard_name'] = 'distance, way of the crow (i.e. direct), to the reference point'
                self.ds['distance'] = da
    def run_Chase2021(self,models_path='../models/'):
    
        """ Method to run Chase et al. (2021) JAMC Neural Network retrieval. THIS NEEDS TENSORFLOW!
        
        Outputs:
        ========
        
        This retrieval will add the following variables to the dataset

        Dml_nn:2d np.array, dim = (along_track,range), unit = mm, name= liquid eq. mass weighted mean diameter 
        Dms_nn: 2d np.array, dim = (along_track,range), unit = mm, name= solid phase mass weighted mean diameter  
        Nw_nn: 2d np.array, dim = (along_track,range), unit = log(m^-4), name= liquid eq. normalized intercept parameter  
        IWC_nn: 2d np.array, dim = (along_track,range), unit = g m^{-3}, name= ice water content 
        
        ========
        """
        import pickle
        import copy

        try:
            import tensorflow as tf
        except:
            print("No Tensorflow. Please go install it with conda")

        with open(models_path+ 'scaler_X.pkl', 'rb') as inp:
            scaler_X = pickle.load(inp)

        with open(models_path+ 'scaler_y.pkl', 'rb') as inp:
            scaler_y = pickle.load(inp)



         

        Ku = copy.deepcopy(self.ds.zFactorFinal[:,:,:,0].values) #use the corrected Ku 
        Ka = copy.deepcopy(self.ds.zFactorMeasured[:,:,:,1].values) #use the raw Ka
        T = copy.deepcopy(self.ds.airTemperature.values) #grab JMA temperature 
        T = T - 273.15 #convert to degC

        #supress warnings. skrews up my progress bar when running in parallel
        def warn(*args, **kwargs):
            pass
        import warnings
        warnings.warn = warn
        #
        
        #load the trained model 
        tf.config.run_functions_eagerly(True)
        model = tf.keras.models.load_model(models_path + 'NN_6by8.h5',custom_objects=None,compile=True)

        #now we have to reshape things to make sure they are in the right shape for the NN model [n_samples,n_features]
        shape_step1 = Ku.shape
        Ku = Ku.reshape([Ku.shape[0]*Ku.shape[1]*Ku.shape[2]])
        Ka = Ka.reshape([Ka.shape[0]*Ka.shape[1]*Ka.shape[2]])
        T = T.reshape([T.shape[0]*T.shape[1]*T.shape[2]])
        
        #make sure that Ku and Ka match 
        Ku[np.isnan(Ka)] = np.nan
        Ka[np.isnan(Ku)] = np.nan
        #
        #mask gates with temperatures > 0 
        Ku[T>0] = np.nan
        #mask gates with DFR < -0.5, shouldnt be ever be lower than this (plus/minus Cal uncert.)
        Ku[(Ku-Ka) < -0.5] = np.nan
        
        #Make sure we only run in on non-nan values to save time 
        ind_masked = np.isnan(Ku) 
        ind_masked2 = np.isnan(Ka)
        Ku_nomask = np.ones(Ku.shape)*-9999.
        Ka_nomask = np.ones(Ka.shape)*-9999.
        T_nomask = np.ones(T.shape)*-9999.
        Ku_nomask[~ind_masked] = Ku[~ind_masked]
        Ka_nomask[~ind_masked] = Ka[~ind_masked]
        T_nomask[~ind_masked] = T[~ind_masked]

        ind = np.where(Ku_nomask!=-9999.)[0]
        
        #scale the input vectors by the mean that it was trained with
        X = np.zeros([Ku_nomask.shape[0],3])
        X[:,0] = (Ku_nomask - scaler_X.mean_[0])/scaler_X.scale_[0] #Ku
        X[:,1] = ((Ku_nomask - Ka_nomask)- scaler_X.mean_[1])/scaler_X.scale_[1] #DFR Ku - Ka
        X[:,2] = (T_nomask - scaler_X.mean_[2])/scaler_X.scale_[2] #T
        #

        #conduct the retrieval 
        yhat = model.predict(X[ind,0:3],batch_size=len(X[ind,0]))
        #scale it properly 
        yhat[:,0] = (yhat[:,0]*scaler_y.scale_[0]) + scaler_y.mean_[0]
        yhat[:,1] = (yhat[:,1]*scaler_y.scale_[1]) + scaler_y.mean_[1]
        yhat[:,2] = (yhat[:,2]*scaler_y.scale_[2]) + scaler_y.mean_[2]
        yhat[:,1] = 10**yhat[:,1] #unlog Dm liquid
        yhat[:,2] = 10**yhat[:,2] #unlog Dm solid

        #fill the Nw array 
        ind = np.where(Ku_nomask!=-9999.)[0]
        Nw = np.zeros(Ku_nomask.shape)
        Nw[ind] = np.squeeze(yhat[:,0])
        #reshape it back to input shape
        Nw = Nw.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
        Nw[Nw==0.0] = np.nan
        #shove it into the xarray dataset 
        self.ds['Nw_nn'] = xr.DataArray(Nw,dims=['nscan','nrayNS','nbin'])

        Dm = np.zeros(Ku_nomask.shape)
        Dm[ind] = np.squeeze(yhat[:,1])
        Dm = Dm.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
        Dm[Dm==0.0] = np.nan 
        self.ds['Dml_nn'] = xr.DataArray(Dm,dims=['nscan','nrayNS','nbin'])
        
        Dm_frozen = np.zeros(Ku_nomask.shape)
        Dm_frozen[ind] = np.squeeze(yhat[:,2])
        Dm_frozen = Dm_frozen.reshape([shape_step1[0],shape_step1[1],shape_step1[2]])
        Dm_frozen[Dm_frozen==0.0] = np.nan
        self.ds['Dms_nn'] = xr.DataArray(Dm_frozen,dims=['nscan','nrayNS','nbin'])

        #calculate IWC 
        Nw = 10**Nw #undo log, should be in m^-4
        Dm = Dm/1000. # convert to m ^4
        IWC = (Nw*(Dm)**4*1000*np.pi)/4**(4) # the 1000 is density of water (kg/m^3)
        IWC = IWC*1000 #convert to g/m^3 
        self.ds['IWC_nn'] = xr.DataArray(IWC,dims=['nscan','nrayNS','nbin'])

        #calculate R (following Chase et al. 2022 paramaterization between log(R) - log(IWC)
        R = 3.64*(IWC**1.06)
        self.ds['R_nn'] = xr.DataArray(R,dims=['nscan','nrayNS','nbin'])
        
        return

