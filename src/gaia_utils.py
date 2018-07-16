"""

Methods to help with the GAIA data.
Python3 compl

HISTORY:
    14.06.2018:
        - first shot
        - adding a class for a dataset from gaia-on-tap
        -
        
"""

__author__  = "SL: ALMA"
__version__ = "0.2.0@2018.07.08"

import math
import numpy as np
import wavelet as wav

DEG2RAD = math.pi / 180.

class gaiaSet:
    
    def __init__(self, data = []):
        "data: GAIA dataset"
        
        self.data = data
        
        
        
    def isHomogeneous(self, tol = 0.1):
        """
        Return True if set is likely homogeneous otherwise False
        The tol is in degrees. projection in RA is checked.
        
        """
        
        
        minDec = min(self.data[ "dec"])
        
        if abs(minDec) > 70:
            longitud = self.data['l']
            latitud  = self.data['b']
        else : 
            longitud = self.data['ra']
            latitud  = self.data['dec'] 
        
        centerlong  = longitud.mean()
        centerlat = latitud.mean()
        
        minlong = min(longitud)
        maxlong = max(longitud)
        minlat = min(latitud)
        maxlat = max(latitud)
        
        halflong  = minlong + (maxlong-minlong)/2.
        halflat = minlat + (maxlat-minlat)/2.
        
        print(halflong)
        print(centerlong)
        print(halflat)
        print(centerlat)
        
        if abs(halflong - centerlong ) < tol  and abs(halflat - centerlat ) < tol:
            return(True)
        else:
            return(False)
        
    def compute_rms(self, image, NSAMPLE = 10):
        "compute the rms using NSAMPLE strip of len(x) / 10"
        
        xdim= image.shape[0]
        ydim= image.shape[1]
                          
        strip = int(xdim / 5)    ## strip length
        
        sample = np.array([])
        
        yrand = np.random.randint(0,ydim,NSAMPLE)
        xrand = np.random.randint(0,xdim-strip-1,NSAMPLE)
        
        for i in range(NSAMPLE):
            testsample = image[xrand[i]:xrand[i]+strip,yrand[i]]
            
            if testsample.std() > 0:
                sample = np.append(sample,testsample)
                
        print(len(sample))
        rms = sample.std()
        return(rms)
            
    
        
    def sampling_filtering(self,x,y,nbin, sigma, LEVELS = 6, xrange = [], yrange = []):
        "Sampling of the points (x,y) and filtering of the images"
        
        
        ### sampling
        if len(xrange) == 0:
            xmin = min(x)
            xmax = max(x)
            ymin = min(y)
            ymax = max(y)
        else:
            xmin = xrange[0]
            xmax = xrange[1]
            ymin = yrange[0]
            ymax = yrange[1]
            
        dx = (xmax-xmin) / nbin
        dy = (ymax-ymin) / nbin
        
        image = np.zeros((nbin,nbin))
        
        ix = ((x - xmin)/dx).astype(int) % nbin
        iy = ((y - ymin)/dy).astype(int) % nbin
        
        image[ix,iy] += 1
        
        
        ## compute the rms by bootstraping
        
        rms = self.compute_rms(image, NSAMPLE = 1000)
        print("## std: %f"%(rms))
        
        ### filtering
        
        wt = wav.wt(verbose = False)
        W = wt.atrous(image,LEVELS)
        Wf = wt.filtering(W,threshold= sigma, waveletNoise = True,imageNoise = rms)
        image_filtered = wt.restore(Wf,0,0)
        
        return(image, image_filtered)
    
    
    
        
        
        
    
