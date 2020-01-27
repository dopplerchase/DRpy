class GPMDPR():

    """
    Author: Randy J. Chase

    This object is intended to help with the efficient processing of GPM-DPR radar files. Currently, xarray cannot
    read the files directly. So here is an attempt to do so. Once in xarray format, the effcient search functions 
    can be used. 
    
    Currently, I do not have this function pass all variables through. Just the files I need to work with.

    """

    def __init__(self,filename=[],do_read =False): 
        self.filename = filename
        self.do_read = do_read
        
    def read(self):
        self.hdf = h5py.File(filename,'r')
        
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

        prh = np.zeros([177,50]) #set up matrix 
        for i in np.arange(0,177): #loop over num range gates
            for j in np.arange(0,50): #loop over scans 
                a = np.arcsin(((re+407)/re)*np.sin(theta[j]))-theta[j] #407 km is the orbit height, re radius of earth, 
                prh[i,j] = (176-(i))*0.125*np.cos(theta[j]+a) #more geometry 
        
        self.height = prh
    
    def toxr(self,snow=True,clutter=True,echotop=True,dropna=True):
    
        """
        
        This function converts the hdf file to an xarray dataset. 
        

        """
        
        lons = self.hdf['NS']['Longitude'][:,12:37]
        lats = self.hdf['NS']['Latitude'][:,12:37]
        
        da = xr.DataArray(self.hdf['MS']['Experimental']['flagSurfaceSnowfall'][:,:], dims=['along_track', 'cross_track'],
                                   coords={'lons': (['along_track','cross_track'],lons),
                                           'lats': (['along_track','cross_track'],lats)})
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
                                   'lats': (['along_track','cross_track'],lats)})
            da.attrs['units'] = 'none'
            da.attrs['standard_name'] = 'flag to remove ground clutter'
            self.xrds['clutter'] = da
        
        da = xr.DataArray(self.hdf['NS']['SLV']['zFactorCorrectedNearSurface'][:,12:37], dims=['along_track', 'cross_track'],
                                   coords={'lons': (['along_track','cross_track'],lons),
                                           'lats': (['along_track','cross_track'],lats)})
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
                                           'lats': (['along_track','cross_track'],lats)})
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'corrected KuPR'
        if clutter:
            da = da.where(self.xrds.clutter==0)
        da = da.where(da >= 12)
        self.xrds['NSKu_c'] = da
        
        da = xr.DataArray(self.hdf['MS']['SLV']['zFactorCorrected'][:,:,:], dims=['along_track', 'cross_track','range'],
                                   coords={'lons': (['along_track','cross_track'],lons),
                                           'lats': (['along_track','cross_track'],lats)})
        da.attrs['units'] = 'dBZ'
        da.attrs['standard_name'] = 'corrected KaPR, MS scan'
        if clutter:
            da = da.where(self.xrds.clutter==0)
        da = da.where(da >= 15)
        self.xrds['MSKa_c'] = da
        
        if echotop:
            self.echotop()
            da = xr.DataArray(self.dummy2, dims=['along_track', 'cross_track','range'],
                                       coords={'lons': (['along_track','cross_track'],self.xrds.lons),
                                               'lats': (['along_track','cross_track'],self.xrds.lats)})
            da.attrs['units'] = 'none'
            da.attrs['standard_name'] = 'flag to remove noise outside cloud/precip top'
            self.xrds['echotop'] = da
        
        da = xr.DataArray(self.hdf['NS']['PRE']['zFactorMeasured'][:,12:37,:], dims=['along_track', 'cross_track','range'],
                                   coords={'lons': (['along_track','cross_track'],lons),
                                           'lats': (['along_track','cross_track'],lats)})
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
                                           'lats': (['along_track','cross_track'],lats)})
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
                                           'lats': (['along_track','cross_track'],lats)})
        da.attrs['units'] = 'mm hr^-1'
        da.attrs['standard_name'] = 'retrieved R, from DPR algo'
        if clutter:
            da = da.where(self.xrds.clutter==0)
        if echotop:
            da = da.where(self.xrds.echotop==0)
        self.xrds['R'] = da
        
        if snow:
            self.xrds = self.xrds.where(self.xrds.flagSurfaceSnow==1,drop=False)
            
         
    def setboxcoords(self,corners=[]):
        
        if len(corners) > 0:
            self.ll_lon = corners[0]
            self.ur_lon = corners[1]
            self.ll_lat = corners[2]
            self.ur_lat = corners[3]
            self.xrds = self.xrds.where((self.xrds.lons >= self.ll_lon) & (self.xrds.lons <= self.ur_lon) & (self.xrds.lats >= self.ll_lat)  & (self.xrds.lats <= self.ur_lat),drop=False)
        else:
            pass