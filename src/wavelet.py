#!/usr/bin/python


"""

Sets of class to provide tools for manipulating Wavelet Transform


HISTORY:

    2010.06.17:
        - first shot to create the main class 
        - Load the shared object for the convolution ("Atrous" algorithm). See the source in IDL/wavelet.
      
      
    2010.06.18:
        - Load convol.so using ctypes
        
    2012.05.29:
        - Change the coding of the A trous WT
        - Add the upscale for the kernel in 1D
        - add filtering (see wavelet.pro for a scaled noise from image to wavelet plane)
        
    2012.05.30:
        - fix the filtering (add noise from wavelet(
        - fix the upscale (scale 1/2)
        
    2012.05.31:
        - Add an argument for NaN to deal with it.
        
        
    2013.12.06:
        - Add the fits class (read, export)
        
    2014.01.02:
        - Add the write fits method
        - still buggy ...
        _ better to look at fitsTools (this part will be deprecated soon)
        
    2014.01.30:
        - fix the filtering method

    2014.01.31:
        - fix the filtering method
        
    2014.05.15:
        - 1D wavelet transform
        
    2014.06.07:
        - change the scipy.convolve to the numpy.convolve to be used in CASA
        
    2017.02.10:
        - comment the io.fits currently.... !!!!!
        
    2017.02.28
        - add verbose option
        
    2018.07.08:
        - update to Python3
        
RUN:
CASA or iPython:
     sys.path.insert(0,'/home/stephane/git/signalanalysis/SignalAnalysis/Wavelet/')
"""


__author__="S. Leon @ ALMA"
__version__="0.4.0@2018.07.08"



import numpy as np
import scipy.signal as signal
import scipy.ndimage as scimage
## import astropy.io.fits as pf



## Kernel  Spline
b3spline1d = np.array([1./16, 1./4, 3./8, 1./4, 1./16])        

class wt:
    """
    Class to perform the WT on images.
    
    """
    
    def __init__(self, verbose = True):
        
        self.verbose = verbose

    def atrous(self,arr, lev, kernel1d=b3spline1d, boundary='symm',NaN= False):
        """
        Do 2d a'trous wavelet transform with B3-spline scaling function
        
        Implementation of the A Trous algorithm.
        - Isotropic and not orthogonal set of wavelet functions (spline)
        - CPU-efficient.
        - See (e.g.) Leon et al., 2000, A&A, 359, 907
        
        
        Entries:
            arr: 2D-image
            lev: Number of planes for the WT
            NaN = False : How to deal with NaN (Default nothing)
            
        """
        
                
        # NaN 
        if NaN :
            image = np.nan_to_num(arr)
        else :
            image = np.copy(arr)
        
        
        
        __x = kernel1d.reshape(1,-1)
        kernel2d = np.dot(__x.T,__x)
        
        if self.verbose:
            print("##WT--A Trous--Plane: %d"%(lev))
        
        approx = signal.convolve2d(image, kernel2d, mode='same', boundary=boundary)  # approximation
        w = image - approx   
                                                                # wavelet details
        if lev <= 0: 
            return image
        
        if lev == 1: 
            if self.verbose:
                print("##WT--A Trous--Plane: 0")
            return [w, approx]
        else:        
            return [w] + self.atrous(approx,lev-1, self.upscale(kernel1d), boundary)
    
    
    def atrous1d(self,arr, lev, kernel1d = b3spline1d,NaN= False):
        """
        Do 1d a'trous wavelet transform with B3-spline scaling function
        
        Implementation of the A Trous algorithm.
        - not orthogonal set of wavelet functions (spline)
        - CPU-efficient.
        - See (e.g.) Leon et al., 2000, A&A, 359, 907
        
        
        Entries:
            arr: 1D-spectrum
            lev: Number of planes for the WT
            NaN = False : How to deal with NaN (Default nothing)
            
        """
        
                
        # NaN 
        if NaN :
            image = np.nan_to_num(arr)
        else :
            image = np.copy(arr)
        
        
        if self.verbose:
            print("##WT--A Trous--Plane: %d"%(lev))
        
        approx = np.convolve(image, kernel1d, mode='same')  
        
        w = image - approx   
        
                                                                
        if lev <= 0: 
            return image
        
        if lev == 1: 
            if self.verbose:
                print("##WT--A Trous--Plane: 0")
            return [w, approx]
        else:        
            return [w] + self.atrous1d(approx,lev-1, self.upscale(kernel1d))
        
        
        
        
    def upscale(self, kernel1d):
        "Upscale (resample) by a factor 2 the kernel.s "
        
        xx = np.arange(len(kernel1d))
        xp = np.arange(0.,len(kernel1d)-0.5,0.5)
        
        nkernel = np.interp(xp,xx,kernel1d)/2.
            
        return(nkernel)
        
        
        
    def filtering(self,wvalue,threshold=3.,mask=(0,0,0,0),waveletNoise = False,imageNoise = 0.):
        """
        Filtering of an image by a clipping of the wavelet coefficients,
        
        Input:
        wvalue : Wavelet Coefficient (list of Wavelet planes)
        threshold : threshold (in sigma) to cut the wavelet coefficients (Default = 3)
        mask : tuple with the box coordinates of the mask (Default = (0,0,0,0))
        waveletNoise : compute the noise of each plane from the image noise
        imageNoise : image STDV (to be  used with waveletNoise)
        """
        
        if self.verbose:
            print("#WT-- Image Filtering")
            print("#WT-- Filtering to be checked")
          
        SIGMA_WAVELET = [0.899677,0.206014,0.0884077,0.0436298,0.0232347,0.0139958,0.00467207]
          
        if mask == (0,0,0,0) and not waveletNoise:
            print("##WT-Filtering--Warning, the mask to compute the noise is (0,0,0,0)")
            
        if waveletNoise and imageNoise == 0.:
            print("##WT-Filtering--Warning, the image noise is 0.")
                  
        wvalueFiltered = []
        nplane = len(wvalue)-1
        indplane = 0
        
        wvalue_c = np.copy(wvalue)
        x1 = int(mask[0])
        y1 = int(mask[2])
        x2 = int(mask[1])
        y2 = int(mask[3])
        
        for plane in wvalue_c:
            planeFiltered = np.copy(plane)
            

            if nplane > 0:
                sigma = np.std(planeFiltered[x1:x2,y1:y2])
                
                if waveletNoise:
                    sigma = imageNoise * SIGMA_WAVELET[indplane]
            
                thresholdPlane = threshold * sigma            
                indT = np.where(abs(planeFiltered) < thresholdPlane)
                         
                if len(indT[0] > 0):
                    planeFiltered[indT[0],indT[1]] = 0.

                if self.verbose:
                    print("##WT--Plane %d Sigma = %e"%(nplane, sigma))
                    print("##WT--Pixel filtered : %d"%(len(indT[0])))
                 
            wvalueFiltered.append(planeFiltered)
            nplane -= 1
            indplane += 1
        
        
        return(wvalueFiltered)
    
            
    def filtering1d(self,wvalue,threshold=3.,mask=(0,0), waveletNoise = False, spectralNoise = 0., sigmaPlane = []):
        """
        Filtering of a spectra by a clipping of the wavelet coefficients,
        
        Input:
        wvalue :       Wavelet Coefficient (list of Wavelet planes)
        threshold :    Threshold (in sigma) to cut the wavelet coefficients (Default = 3)
        mask :         Tuple with the box coordinates of the mask (Default = (0,0,0,0))
        waveletNoise : Compute the noise of each plane from the spectra noise
        spectraNoise :   spectra STDV (to be  used with waveletNoise)
        sigmaPlane:    provide the array of the sigma for each plane
        """
        
        if self.verbose:
            print("#WT--Spectrum Filtering")
          
            
        SIGMA_WAVELET = [0.899677,0.206014,0.0884077,0.0436298,0.0232347,0.0139958,0.00467207]
          
        if mask == (0,0) and not waveletNoise:
            print("##WT-Filtering--Warning, the mask to compute the noise is (0,0)")
            
        if waveletNoise and spectralNoise == 0.:
            print("##WT-Filtering--Warning, the image noise is 0.")
                  
        wvalueFiltered = []
        nplane = len(wvalue)-1
        indplane = 0
        
        wvalue_c = np.copy(wvalue)
        x1 = int(mask[0])
        x2 = int(mask[1])
        
        sigmaProvided = False
        
        if len(sigmaPlane) > 0:
            sigmaProvided = True
            sigmaPlane.reverse()
        
        for plane in wvalue_c:
            planeFiltered = np.copy(plane)
            
            if nplane > 0:
                
                if sigmaProvided:
                    sigma = sigmaPlane[nplane-1]
                elif mask != (0,0) :
                    sigma = np.std(planeFiltered[x1:x2])
                
                if waveletNoise:
                    sigma = spectralNoise * SIGMA_WAVELET[indplane]
            
                thresholdPlane = threshold * sigma            
                indT = np.where(abs(planeFiltered) < thresholdPlane)
                    
                if len(indT[0] > 0):
                    planeFiltered[indT[0]] = 0.

                if self.verbose:
                    print("##WT--Plane %d Sigma = %e"%(nplane, sigma))
                    print("##WT--Pixel filtered : %d"%(len(indT[0])))
                 
            wvalueFiltered.append(planeFiltered)
            nplane -= 1
            indplane += 1
        
        
        return(wvalueFiltered)       
    def restore(self,wvalue,plane1,plane2):
        
        """
          Return an image by summing plane1 to plane2 (if plane1 = plane2 = 0 sum all the planes)
        """      
        
        if self.verbose:
            print("##WT--Restore-plane: %d to %d"%(plane1,plane2))
        
        if plane1 == 0 and plane2 ==0:
            plane1 = 0
            plane2 = len(wvalue)-1
        
        image = np.copy(wvalue[plane1])
        
        image[:,:] = 0.
          
        for i in range(plane1,plane2+1):
            if self.verbose:
                print("##WT--Restore-plane: %d"%(i))
            image += wvalue[i]
        
        return(image)
                 
       
    def restore1d(self,wvalue,plane1,plane2):
        
        """
          Return an array  by summing plane1 to plane2 (if plane1 = plane2 = 0 sum all the planes)
        """      
        
        if self.verbose:
            print("##WT--Restore-plane: %d to %d"%(plane1,plane2))
        
        if plane1 == 0 and plane2 ==0:
            plane1 = 0
            plane2 = len(wvalue)-1
        
        image = np.copy(wvalue[plane1])
        
        image[:] = 0.
          
        for i in range(plane1,plane2+1):
            if self.verbose:
                print("##WT--Restore-plane: %d"%(i))
            image += wvalue[i]
        
        return(image)             
                
class fits:
    """ class to read and write fits file to be treated by WT.
    """
    
    
    def __init__(self):
        
        print("# Deprecated by fitsTools.")
    
    
    def read(self,fitsname):
        """ Read a fitsname. The default is set to a standard ALMA  image """
        
        
        hdulist = pf.open(fitsname)        
        arr = hdulist[0].data[0][0]
        hdulist.close()
        
        print("##FITS--read %s"%(fitsname))
        
        return(arr, hdulist)
    
    
    
    def write(self,fitsname, arr, hd = None):
        "Write a fits file. If hd is none, it does not write any header"
        
        if hd !=  None:
            pf.writeto(fitsname, arr, header = hd[0].header)
        else:
            pf.writeto(fitsname, arr)
        
        print("##FITS--write %s"%(fitsname))
        print("##FITS--write WCS not recognized by CASA but by ds9...")
        
    
        
        
        
    
    
    
        
    
        
