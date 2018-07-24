"""

Methods to help with the GAIA data.
Python3 compl

HISTORY:
    14.06.2018:
        - first shot
        - adding a class for a dataset from gaia-on-tap
    
    19.07.2018:
        - adding filterinf for data
        - gaia query
        - color and magnitude
        
    24.07.2018
        - changing the normalization 
        
"""

__author__  = "SL, QV: ALMA"
__version__ = "0.3.1@2018.07.24"

# Suppress warnings
import warnings
warnings.filterwarnings('ignore')

import math , shutil
import numpy as np
import wavelet as wav

import astropy.coordinates as coord

from astropy.io.votable import parse
from astropy.table import Table
from astropy import units as u

from astroquery.gaia import Gaia

DEG2RAD = math.pi / 180.


DIMMAX = 10

class source:
    "Class to download gaia data for a source"
    
    
    ########################
    def __init__(self, name):
        
        self.name = name
        self.weight = np.ones(DIMMAX)
    
    
    ################################
    def query(self, radius, errtol = 0.1, dump = False,table="gaiadr2.gaia_source"):
        "do a conesearch"
        
        c = coord.SkyCoord.from_name(self.name)
        c_icrs = coord.SkyCoord(ra= c.ra.deg * u.degree, dec= c.dec.deg  *u.degree, frame='icrs')

        pos = c_icrs.galactic
        self.l_cluster = pos.l.value 
        self.b_cluster = pos.b.value
        self.ra  = c.ra.deg
        self.dec = c.dec.deg
        
        
        queryaql = "SELECT * FROM {table} WHERE CONTAINS(POINT('ICRS',{table}.ra,{table}.dec),  \
                                    CIRCLE('ICRS',{ra:.10f},{dec:.10f},{radius:.10f})) = 1  AND abs(pmra_error/pmra)<{err:10f}  AND abs(pmdec_error/pmdec)< {err:.10f} AND abs(parallax_error/parallax)< {err:.10f};".format(table=table, ra=c.ra.deg, dec=c.dec.deg,radius=radius, err=errtol)
        
        print(queryaql)
        job = Gaia.launch_job_async(queryaql, dump_to_file=dump)
        self.data = job.get_results()

        if dump:
            filename = job.get_output_file()
            filedst = "%s-%3.1fdeg.vot"%(self.name, radius)
            shutil.move(filename,filedst)
            print("## %s created"%(filedst))
        else:
            filename = None
            
        print("## Query for %s done"%(self.name))
        print("## Total stars: %d"%(len(self.data)))
        
        return(filename)
    
    
    ########################
    def read_votable(self, voname):
        "rad a votable"
        

        votable = parse(voname)
        for table in votable.iter_tables():
            data = table.array
        #print(data.dtype.names)
    
        self.data = data
        
        print("## %s read..."%(voname))
        print("## Total stars: %d"%(len(data)))
              
        return(len(data))
    
    ###################################################################################################
    def convert_filter_data(self, dist_range = [0., 2000], vra_range = [-200,200], vdec_range = [-200.,200]) :
        
        lgal = self.data['l']
        bgal = self.data['b']
        pmas = self.data['parallax']
        distance = 1000. / np.ma.filled(pmas, -999999.)    #distance = 1000/parallax
        pmra = np.ma.filled(self.data['pmra'], -9999999.)   # PM RA
        pmdec= np.ma.filled(self.data['pmdec'],-9999999.)   # PM Dec
        vdec = 4.74 * pmdec / pmas   ##? (pour )
        vra  = 4.74 * pmra  / pmas   # pour avoir des km.s-1
    
        
    
        #para_abs_error = data['parallax_error']/data['parallax']
        g  =  np.ma.filled(self.data['phot_g_mean_mag'], 99.)
        bp =  np.ma.filled(self.data['phot_bp_mean_mag'], 999.)
        rp =  np.ma.filled(self.data['phot_rp_mean_mag'], 9999.)
    
        #filtering
    
        i1 = np.where((distance >= dist_range[0]) & (distance < dist_range[1]))
        i2 = np.where((vra >= vra_range[0]) & (vra < vra_range[1]))
        i12 = np.intersect1d(i1,i2)
        i3 = np.where((vdec >= vdec_range[0]) & (vdec < vdec_range[1]))
        ifinal = np.intersect1d(i12,i3)
        
        gbar = g[ifinal] + 5*np.log10(pmas[ifinal]) + 2
        
        self.df = np.array([lgal[ifinal],bgal[ifinal], distance[ifinal], vra[ifinal], vdec[ifinal], gbar, g[ifinal]-rp[ifinal], bp[ifinal]-g[ifinal]]).T
        
        print("## Conversion done...")
        print("## Stars selected: %d"%(len(ifinal)))
        return()
    
    
    ###############################
    def normalization_normal(self):
        "normalize to the STD and MEAN"
        
        self.dfnorm = np.zeros(self.df.shape)
    
        for i in range(self.df.shape[1]) :
            
            self.dfnorm[:,i] = self.weight[i] * ( self.df[:,i] - np.mean(self.df[:,i])) / np.std(self.df[:,i]) 
        
        print("## Normalization done on filtered data..")        
        return()
 
 
 
    ######################
    #Normalized the 5d datask with linear projection from [min,max] to [0,1]
    def normalization_minmax(self):
        
        self.dfnorm = np.zeros(self.df.shape)
        self.normalization_vector = np.zeros((DIMMAX,2)) #Represente max and min    
    
        for i in range(self.df.shape[1]) :
            self.normalization_vector[i,0] = np.max(self.df[:,i]) # max
            self.normalization_vector[i,1] = np.min(self.df[:,i]) # min
            self.dfnorm[:,i] = self.weight[i]*(self.df[:,i]-self.normalization_vector[i,1])/(self.normalization_vector[i,0]-self.normalization_vector[i,1])  
        
        print("## Normalization done on filtered data..")
        return()

    ###############################################  
    def convert_to_cartesian(self, offCenter = []):
        "Convert ra,dec (ICRS) and distance (pc) to Cartesian reference. Off is the offset in Lgal,Bgal"
    
        xx = np.zeros(len(self.df[:,0]))
        yy = np.zeros(len(self.df[:,0]))
        zz = np.zeros(len(self.df[:,0]))
    
        if len(offCenter) == 0:
            offCenter[0] = self.l_cluster
            offCenter[1] = self.b_cluster
        
        lgalOff = self.df[:,0] - offCenter[0]
        bgalOff = self.df[:,1] - offCenter[1]
    
    
        for i in range(len(lgalOff)):
            c = coord.SkyCoord(l=lgalOff[i]*u.degree, b=bgalOff[i]*u.degree, distance=self.df[:,2] *u.pc, frame='galactic')
        
            xx[i] = c.cartesian.x.value
            yy[i] = c.cartesian.y.value
            zz[i] = c.cartesian.z.value
        
    
        return(xx, yy, zz)


##################################################################################################################
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
    
    
    
        
        
        
    
