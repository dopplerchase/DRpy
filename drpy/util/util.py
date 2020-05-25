from __future__ import absolute_import

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def percentile(n):
    def percentile_(x):
        return np.percentile(x, n)
    percentile_.__name__ = 'percentile_%s' % n
    return percentile_

def boxbin(x,y,xedge,yedge,c=None,figsize=(5,5),cmap='viridis',mincnt=10,vmin=None,vmax=None,edgecolor=None,powernorm=False,
           ax=None,normed=False,method='mean',quantile=None,alpha=1.0,cbar=True):
    
    """ This function will grid data for you and provide the counts if no variable c is given, or the median if 
    a variable c is given. In the future I will add functionallity to do the median, and possibly quantiles. 
    
    x: 1-D array 
    y: 1-D array 
    xedge: 1-D array for xbins 
    yedge: 1-D array for ybins
    
    c: 1-D array, same len as x and y 
    
    returns
    
    axis handle 
    cbar handle 
    C matrix (counts or median values in bin)
    
    """
    
    midpoints = np.empty(xedge.shape[0]-1)
    for i in np.arange(1,xedge.shape[0]):
        midpoints[i-1] = xedge[i-1] + (np.abs(xedge[i] - xedge[i-1]))/2.
    
    #note on digitize. bin 0 is outside to the left of the bins, bin -1 is outside to the right
    ind1 = np.digitize(x,bins = xedge) #inds of x in each bin
    ind2 = np.digitize(y,bins = yedge) #inds of y in each bin

    

    if c is None:
        c = np.zeros(len(x))
        df = pd.DataFrame({'x':ind1,'y':ind2,'c':c})
        df2 = df.groupby(["x","y"]).count()
        df = df2.where(df2.values >= mincnt).dropna()
        C = np.ones([xedge.shape[0]-1,yedge.shape[0]-1])*-9999
        for i,ii in enumerate(df.index.values):
            #dont plot points outside the bounds
            if (ii[0]==0) or (ii[1] == 0) or (ii[0] >= len(xedge)-1) or (ii[1] >= len(yedge)-1):
                pass
            else:
                #shift indexs so that bin 1 is actually bin 0. 
                C[ii[0]-1,ii[1]-1] = df.c.values[i]
        C = np.ma.masked_where(C == -9999,C)
        
        if normed:
            n_samples = np.ma.sum(C)
            C = C/n_samples
            C = C*100
            print('n_samples= {}'.format(n_samples))
        
        if ax is None:
            fig = plt.figure(figsize=(5,5))
            ax = plt.gca()
        else:
            pass
            
        if powernorm:
            pm = ax.pcolormesh(xedge,yedge,C.transpose(),cmap=cmap,edgecolor=edgecolor,norm=colors.PowerNorm(gamma=0.5),vmin=vmin,vmax=vmax,alpha=alpha)
            
            if cbar:
                cbar = plt.colorbar(pm,ax=ax)
        else:
            pm = ax.pcolormesh(xedge,yedge,C.transpose(),cmap=cmap,vmin=vmin,vmax=vmax,edgecolor=edgecolor,alpha=alpha)
            if cbar:
                cbar = plt.colorbar(pm,ax=ax)
            
        return ax,cbar,C
    
    else:
        df = pd.DataFrame({'x':ind1,'y':ind2,'c':c})
        if method=='mean':
            df2 = df.groupby(["x","y"])['c'].mean()
        elif method=='std':
            df2 = df.groupby(["x","y"])['c'].std()
        elif method=='median':
            df2 = df.groupby(["x","y"])['c'].median()
        elif method=='qunatile':
            if quantile is None:
                print('No quantile given, defaulting to median')
                quantile = 0.5
            else:
                pass
            df2 = df.groupby(["x","y"])['c'].apply(percentile(quantile*100))
            
            
        df3 = df.groupby(["x","y"]).count()
        df2 = df2.to_frame()
        df2.insert(1,'Count',df3.values)
        df = df2.where(df2.Count >= mincnt).dropna()
        C = np.ones([xedge.shape[0]-1,yedge.shape[0]-1])*-9999
        for i,ii in enumerate(df.index.values):
            #dont plot points outside the bounds
            if (ii[0]==0) or (ii[1] == 0) or (ii[0] >= len(xedge)-1) or (ii[1] >= len(yedge)-1):
                pass
            else:
                #shift indexs so that bin 1 is actually bin 0.
                C[ii[0]-1,ii[1]-1] = df.c.values[i]

        C = np.ma.masked_where(C == -9999,C)

        if ax is None:
            fig = plt.figure(figsize=(5,5))
            ax = plt.gca()
        else:
            pass
        
        if powernorm:
            pm = ax.pcolor(xedge,yedge,C.transpose(),cmap=cmap,vmin=vmin,vmax=vmax,norm=colors.PowerNorm(gamma=0.5),alpha=alpha)
            if cbar:
                cbar = plt.colorbar(pm,ax=ax)
        else:
            
            pm = ax.pcolor(xedge,yedge,C.transpose(),cmap=cmap,vmin=vmin,vmax=vmax,alpha=alpha)
            if cbar: 
                cbar = plt.colorbar(pm,ax=ax)
            
    return ax,cbar,C
