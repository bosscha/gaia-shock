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
        
    07.08.21018:
        - adding coor option to query
        
    11.08.2018:
        - add a method to extract DBSCAN clusters.
        
    14.08.2018:
        - adding cartesian conversion to df and updating names
        
    15.08.2018:
        - fixing bug
        - adding normalizationVecto to make the definition more general
        
    17.08.2918;
        - fixing bugs
        
    18.08.2018:
        - minmax normalisation in perBlocks
        - normalization per block scaled down to sum of weight
        - update dbscan_ results for cartesian
        

"""

__author__  = "SL, QV: ALMA"
__version__ = "0.4.2@2018.08.18"

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

from sklearn import cluster

DEG2RAD = math.pi / 180.


DIMMAX = 10

class source:
    "Class to download gaia data for a source and to perform a clustering extractopn"
    
    
    ########################
    def __init__(self, name):
        
        self.name       = name
        self.weight     = np.ones(DIMMAX)
        self.weightcart = np.ones(DIMMAX)
        self.cartesuan  = False
    
    
    ################################
    def query(self, radius, coordCluster = [], errtol = 0.2, dump = False,table="gaiadr2.gaia_source"):
        "do a conesearch"
        
        self.weighted = False
        
        if len(coordCluster) == 2:
            c = coord.SkyCoord(ra=coordCluster[0]*u.degree, dec=coordCluster[1]*u.degree, frame='icrs')
            c_icrs = c
        else:
            try:
                c = coord.SkyCoord.from_name(self.name)
                c_icrs = coord.SkyCoord(ra= c.ra.deg * u.degree, dec= c.dec.deg  *u.degree, frame='icrs')
            except:
                print("## Cluster name not found ...")
            
            
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
        
        return(filedst)
    
    
    ########################
    def read_votable(self, voname):
        "rad a votable"
        
        self.weighted = False
        
        votable = parse(voname)
        for table in votable.iter_tables():
            data = table.array
        #print(data.dtype.names)
    
        self.data = data
        
        print("## %s read..."%(voname))
        print("## Total stars: %d"%(len(data)))
              
        return(len(data))
    
    ###################################################################################################
    def filter_data(self, dist_range = [0., 2000], vra_range = [-200,200], vdec_range = [-200.,200], mag_range =[-1e9, 1e9]) :
        
        self.weighted = False
        
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
        i4 = np.intersect1d(i12,i3)
        i5 = np.where((g >= mag_range[0]) & (g < mag_range[1]) & (bp >= mag_range[0]) & (bp < mag_range[1]) & (rp >= mag_range[0]) & (rp < mag_range[1]))
        ifinal = np.intersect1d(i4,i5)
        
        
        gbar = g[ifinal] + 5*np.log10(pmas[ifinal]) + 2
        
        self.df = np.array([lgal[ifinal],bgal[ifinal], distance[ifinal], vra[ifinal], vdec[ifinal], gbar, g[ifinal]-rp[ifinal], bp[ifinal]-g[ifinal]]).T
        
        print("## Conversion done...")
        print("## Stars selected: %d"%(len(ifinal)))
        return()
    
    
    
    ###############################
    def normalization_normal(self):
        "normalize to the STD and MEAN"
        
        self.dfnorm = np.zeros(self.df.shape)
        if self.cartesian:
            self.dfcartnorm = np.zeros(self.df.shape)
        self.weighted = True
        
        for i in range(self.df.shape[1]) :
            self.dfnorm[:,i] = self.weight[i] * ( self.df[:,i] - np.mean(self.df[:,i])) / np.std(self.df[:,i]) 
            if self.cartesian:
               self.dfcartnorm[:,i] = self.weightcart[i] * ( self.df[:,i] - np.mean(self.dfcart[:,i])) / np.std(self.dfcart[:,i]) 
        
        print("## Normalization Normal-Gauss done on filtered data..")        
        return()
    
 
    ################################
    def  normalization_PerBlock(self, block, weightblock, cartesian = False, norm = "Identity", density = False):
        """
        To apply the same weight on subset (i axis). Typically Spatial, velocity and magnitudes
        blocks is a list of index list to be gathered
        If cartesian = True only for cartesian (should be added before)
        If density = True normalized as well by number of stars
        """
        
        totalWeight = np.sum(weightblock)
        
        if not cartesian:
            self.dfnorm = np.zeros(self.df.shape)
        elif cartesian:
            self.dfcartnorm = np.zeros(self.df.shape)
        self.weighted = True  
        
        for axis, weight in zip(block, weightblock):
            if not cartesian:
                normK = self.normalizationVector(norm, density, self.df[:,axis])
                normK[1] = normK[1] / totalWeight
                self.dfnorm[:,axis]    =   weight * ( self.df[:,axis] - normK[0] ) / normK[1]
            elif cartesian:
                normK = self.normalizationVector(norm, density, self.dfcart[:,axis])  
                normK[1] = normK[1] / totalWeight               
                self.dfcartnorm[:,axis] =  weight * ( self.dfcart[:,axis]  - normK[0] ) / normK[1] 
                
            
        print("## Normalization (weighted) per block done on %s filtered data..."%("cartesian" if cartesian else ""))
        return()
    
    
    #########################################
    def normalizationVector(self, norm, density, arr):
        "return a normalisation vector"
        
        vecNorm = [0.0, 1.0]
        
        if norm == "identity":
            vecNorm = [0.0, 1.0]
            
        if norm == "averagestep":
            sortArr = np.sort(arr, axis = None)
            diffmean = np.median(np.diff(sortArr))
            vecNorm  = [ 0.0, diffmean]
            
        if norm == "normal":
            stdArr = np.std(arr, axis = None)
            meanArr = np.mean(arr, axis = None)
            vecNorm  = [meanArr , stdArr]
            
        if norm == "minmax":
            minarr  = np.min(arr)
            maxarr  = np.max(arr)
            vecNorm = [minarr, maxarr-minarr]
            
        if density:
            vecNorm [1] = vecNorm[1] * len(arr[:,0])
            
        return(vecNorm)
        
 
    ######################
    #Normalized the 5d datask with linear projection from [min,max] to [0,1]
    def normalization_minmax(self):
        
        self.dfnorm = np.zeros(self.df.shape)
        self.normalization_vector = np.zeros((DIMMAX,2)) #Represente max and min    
        self.weighted = True
    
        for i in range(self.df.shape[1]) :
            self.normalization_vector[i,0] = np.max(self.df[:,i]) # max
            self.normalization_vector[i,1] = np.min(self.df[:,i]) # min
            self.dfnorm[:,i] = self.weight[i]*(self.df[:,i]-self.normalization_vector[i,1])/(self.normalization_vector[i,0]-self.normalization_vector[i,1])  
        
        print("## Normalization minmax done on filtered data..")
        
        
        return()

    ###############################################  
    def add_cartesian(self, offCenter = []):
        "Create a dfcart by converting l,b (ICRS) and distance (pc) to Cartesian reference. Off is the offset in Lgal,Bgal"
    
        self.dfcart = np.copy(self.df)
    
        if len(offCenter) == 0:
            offCenter = [0.,0.]
            offCenter[0] = np.mean(self.df[:,0])
            offCenter[1] = np.mean(self.df[:,1])
        
        lgalOff = self.df[:,0] - offCenter[0]
        bgalOff = self.df[:,1] - offCenter[1]
    
    
        coc = coord.SkyCoord(l=lgalOff*u.degree, b=bgalOff*u.degree, distance=self.df[:,2] *u.pc, frame='galactic')
            
        self.dfcart[:,0] = coc.cartesian.x.value
        self.dfcart[:,1] = coc.cartesian.y.value
        self.dfcart[:,2] = coc.cartesian.z.value
        
        self.cartesian = True
    
        print("## Cartesian coordinates added...")
        return(True)


    #######################################################
    def dbscan_(self, eps, min_samples, cartesian = False):
        "Perform the dbscan on dfnorm. If not weighted returns error. If cartesian uses dfcartnorm"
        
        
        if not self.weighted:
            print("## DBSCAN error, data not weighted..")
            return([])
        
        dbscan = cluster.DBSCAN(eps = eps, min_samples = min_samples)
        if cartesian:
            dbscan.fit(self.dfcartnorm)
        else:
            dbscan.fit(self.dfnorm)
            
        labels_ = dbscan.labels_
        unique_labels = set(labels_)
        n_clusters_ = len(set(labels_)) - (1 if -1 in labels_ else 0)
        
        if cartesian:
            result_ = {}
            result_['label'] = []
            result_['nstars'] = []  
            result_['pos'] = []
            result_['pos_std'] = []
            result_['vel'] = []
            result_['vel_std'] = []
        else:
            result_ = {}
            result_['label'] = []
            result_['nstars'] = []
            result_['distance'] = []
            result_['distance_std'] = []    
            result_['pos'] = []
            result_['pos_std'] = []
            result_['vel'] = []
            result_['vel_std'] = []  
    
        for i in range(-1,n_clusters_):
            if cartesian:
                result_['label'].append(i)
                result_['nstars'].append(len(labels_[np.where(labels_ == i)]))
                result_['pos'].append([np.median(self.dfcart[np.where(labels_ == i),0]), np.median(self.dfcart[np.where(labels_ == i),1]),np.median(self.dfcart[np.where(labels_ == i),2])])
                result_['pos_std'].append([np.std(self.dfcart[np.where(labels_ == i),0]), np.std(self.dfcart[np.where(labels_ == i),1]),np.std(self.dfcart[np.where(labels_ == i),2])]) 
                result_['vel'].append([np.median(self.dfcart[np.where(labels_ == i),3]), np.median(self.dfcart[np.where(labels_ == i),4])])
                result_['vel_std'].append([np.std(self.dfcart[np.where(labels_ == i),3]), np.std(self.dfcart[np.where(labels_ == i),4])])         

            else:
                result_['label'].append(i)
                result_['nstars'].append(len(labels_[np.where(labels_ == i)]))
                result_['distance'].append(np.median(self.df[np.where(labels_ == i),2]))
                result_['distance_std'].append(np.std(self.df[np.where(labels_ == i),2]))
                result_['pos'].append([np.median(self.df[np.where(labels_ == i),0]), np.median(self.df[np.where(labels_ == i),1])])
                result_['pos_std'].append([np.std(self.df[np.where(labels_ == i),0]), np.std(self.df[np.where(labels_ == i),1])])
                result_['vel'].append([np.median(self.df[np.where(labels_ == i),3]), np.median(self.df[np.where(labels_ == i),4])])
                result_['vel_std'].append([np.std(self.df[np.where(labels_ == i),3]), np.std(self.df[np.where(labels_ == i),4])])         

    
        return(labels_, result_) 
        
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
    
    
    
        
        
        
    
