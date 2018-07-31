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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import cluster

DEG2RAD = math.pi / 180.


DIMMAX = 8

class source:
    "Class to download gaia data for a source"
    
    
    ########################
    def __init__(self, name, radius, errtol):

        self.name = name        
        self.radius = radius
        self.errtol = errtol
        self.weight = np.ones(DIMMAX)
        self.data_name = ['distance','lgal','bgal','vdec','vra',r'$G + 5 * log_{10}\bar{\omega} + 2$','$G - R_p$','$B_p - G$']
        self.data_name_cart = ['distance (x)','y','z','vdec','vra',r'$G + 5 * log_{10}\bar{\omega} + 2$','$G - R_p$','$B_p - G$']
    
    
    ################################
    def set_weight(self, weight) :
        "Set the variable weight and scale the vector with its highest value"
        
        self.weight = weight/np.max(weight)    
    
    ################################
    def query(self, dump = False,table="gaiadr2.gaia_source"):
        "do a conesearch"
        
        c = coord.SkyCoord.from_name(self.name)
        c_icrs = coord.SkyCoord(ra= c.ra.deg * u.degree, dec= c.dec.deg  *u.degree, frame='icrs')

        pos = c_icrs.galactic
        self.l_cluster = pos.l.value 
        self.b_cluster = pos.b.value
        self.ra  = c.ra.deg
        self.dec = c.dec.deg
        
        
        queryaql = "SELECT * FROM {table} WHERE CONTAINS(POINT('ICRS',{table}.ra,{table}.dec),  \
                                    CIRCLE('ICRS',{ra:.10f},{dec:.10f},{radius:.10f})) = 1  \
                                    AND abs(pmra_error/pmra)<{err:10f}  AND abs(pmdec_error/pmdec)< {err:.10f} \
                                    AND abs(parallax_error/parallax)< {err:.10f};".format(table=table, ra=c.ra.deg, dec=c.dec.deg,radius=self.radius, err=self.errtol)
        
        print(queryaql)
        job = Gaia.launch_job_async(queryaql, dump_to_file=dump)
        self.data = job.get_results()

        if dump:
            filename = job.get_output_file()
            filedst = "%s-%3.1fdeg-%serr.vot"%(self.name, self.radius, self.errtol)
            shutil.move(filename,filedst)
            print("## %s created"%(filedst))
        else:
            filename = None
            
        print("## Query for %s done"%(self.name))
        print("## Total stars: %d"%(len(self.data)))
        
        return(filename)
    
    
    ########################
    def read_votable(self, voname = None):
        "rad a votable"
        
        # if no voname, find it with name, raduis, errtol
        if voname == None :
            voname = "%s-%3.1fdeg-%serr.vot"%(self.name, self.radius, self.errtol)
            
        
        print("## %s read..."%(voname))

        votable = parse(voname)
        for table in votable.iter_tables():
            data = table.array
        #print(data.dtype.names)
    
        self.data = data
        
        print("## Total stars: %d"%(len(data)))
              
        return(len(data))
    
    ###################################################################################################
    def convert_filter_data(self, dist_range = [0., 2000], vra_range = [-200,200], vdec_range = [-200.,200], without_mag = True) :
        
        lgal = self.data['l']
        bgal = self.data['b']
        pmas = self.data['parallax']
        distance = 1000. / np.ma.filled(pmas, -999999.)    #distance = 1000/parallax
        pmra = np.ma.filled(self.data['pmra'], -9999999.)   # PM RA
        pmdec= np.ma.filled(self.data['pmdec'],-9999999.)   # PM Dec
        vdec = 4.74 * pmdec / pmas   ##? (pour )
        vra  = 4.74 * pmra  / pmas   # pour avoir des km.s-1
    
        
        g  =  np.ma.filled(self.data['phot_g_mean_mag'], 99.)
        bp =  np.ma.filled(self.data['phot_bp_mean_mag'], 999.)
        rp =  np.ma.filled(self.data['phot_rp_mean_mag'], 9999.)
    
        #filtering
    
        i1 = np.where((distance >= dist_range[0]) & (distance < dist_range[1]))
        i2 = np.where((vra >= vra_range[0]) & (vra < vra_range[1]))
        i12 = np.intersect1d(i1,i2)
        i3 = np.where((vdec >= vdec_range[0]) & (vdec < vdec_range[1]))
        ifinal = np.intersect1d(i12,i3)
        
        # if we don't consider masked magnitude
        if without_mag :
            i1 = np.where(g < 99.)[0]
            i2 = np.where(bp < 99.)[0]
            i12 = np.intersect1d(i1,i2)
            i3 = np.where(rp < 99.)[0]
            i4 = np.intersect1d(i12,i3)
            ifinal = np.intersect1d(i4,ifinal)
            
        print()
        
        gbar = g[ifinal] + 5*np.log10(pmas[ifinal]) + 2
        
        self.df = np.array([distance[ifinal], lgal[ifinal], bgal[ifinal], vra[ifinal], vdec[ifinal], gbar, g[ifinal]-rp[ifinal], bp[ifinal]-g[ifinal]]).T
        
        print("## Conversion done...")
        print("## Stars selected: %d"%(len(ifinal)))
    
    
    ###############################
    def normalization_normal(self):
        "normalize to the STD and MEAN"
        
        self.dfnorm = np.zeros(self.df.shape)
    
        for i in range(self.df.shape[1]) :
            
            self.dfnorm[:,i] = self.weight[i] * ( self.df[:,i] - np.mean(self.df[:,i])) / np.std(self.df[:,i]) 
        
        print("## Normalization done on filtered data..")
 
 
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

    ###############################################  
    def convert_to_cartesian(self, centering = True):
        #conversion of distance, lgal, bgal into cartesian coordinate
        #centering True = x coordinate point to cluster center
        "Convert ra,dec (ICRS) and distance (pc) to Cartesian reference. Off is the offset in Lgal,Bgal"
        
        xx = np.zeros(self.df.shape[0]);  yy = np.zeros(self.df.shape[0]);  zz = np.zeros(self.df.shape[0])
        
        dist = np.copy(self.df[:,0]);    lgal = np.copy(self.df[:,1]);    bgal = np.copy(self.df[:,2])
        
        if centering :
            lgal = lgal - np.mean(lgal)
            bgal = bgal - np.mean(bgal)
        
        
        for i in range(len(lgal)):
            c = coord.SkyCoord(l=lgal[i]*u.degree, b=bgal[i]*u.degree, distance=dist[i]*u.pc, frame='galactic')
            
            xx[i] = c.cartesian.x.value
            yy[i] = c.cartesian.y.value
            zz[i] = c.cartesian.z.value
            
        self.dfcart = np.copy(self.df)
        self.dfcart[:,:3] = np.array([xx,yy,zz]).T


    ##############################################
    def HRD(self, size=0.1, colorbar = True):
        "Plot the HD diagram"
        
        plt.figure(figsize=(15,6))
        G_max = np.max(self.df[:,5])
        G_min = np.min(self.df[:,5])
        
        for i in (6,7) :
            plt.subplot(1,2,i-5)
            plt.ylim(G_max, G_min)
            plt.title(self.name, fontsize=20)
            if colorbar : 
                plt.scatter(self.df[:,i], self.df[:,5], s=size, c=self.df[:,0], cmap='gist_stern')
                clb = plt.colorbar()
                clb.set_label('distance', labelpad=-40, y=1.05, rotation=0)
            else : plt.scatter(self.df[:,i], self.df[:,5], s=size, c='k')
            plt.xlabel(self.data_name[i], fontsize=18)
            if i == 6 : plt.ylabel(self.data_name[5], fontsize=22)
        
        plt.show()

    ##############################################
    def plot_information(self, size=0.1, cartesian=False, HRD=True, ilabel=[]) :
        "Plot some graphs about data"
        
        plt.figure(figsize=(19,19))                
        for i_x, i_y, i in zip((1,2,1,1,1,3),(0,0,2,3,4,4),(1,2,3,4,5,6)) :
            plt.subplot(3,2,i)
            if cartesian : 
                plt.scatter(self.dfcart[:,i_x],self.dfcart[:,i_y],s=size,c='k')
                plt.scatter(self.dfcart[ilabel,i_x],self.dfcart[ilabel,i_y],s=size*20,c='r')  
                plt.xlabel(self.data_name_cart[i_x], fontsize=25)
                plt.ylabel(self.data_name_cart[i_y], fontsize=25)              
            else :
                plt.scatter(self.df[:,i_x],self.df[:,i_y],s=size,c='k')
                plt.scatter(self.df[ilabel,i_x],self.df[ilabel,i_y],s=size*20,c='r')
                plt.xlabel(self.data_name[i_x], fontsize=25)
                plt.ylabel(self.data_name[i_y], fontsize=25)
        plt.show()
        if HRD : self.HRD(size)
        
    ##############################################
    def plot_3D(self, size=0.1, cartesian=False, axes = (0,1,2), ilabel=[]) :
        "Plot in 3D the 3 axis"

        fig = plt.figure(figsize=(15,10))
        ax = fig.add_subplot(111, projection='3d')
        if cartesian == False : 
            ax.scatter(self.df[:,axes[0]], self.df[:,axes[1]], self.df[:,axes[2]], zdir='z', s=size, c='k', depthshade=True)
            ax.scatter(self.df[ilabel,axes[0]], self.df[ilabel,axes[1]], self.df[ilabel,axes[2]], zdir='z', s=size*30, c='r', depthshade=True)
            ax.set_xlabel(self.data_name[axes[0]], fontsize=35); ax.set_ylabel(self.data_name[axes[1]], fontsize=25); ax.set_zlabel(self.data_name[axes[2]], fontsize=25)
        else :
            ax.scatter(self.dfcart[:,axes[0]], self.dfcart[:,axes[1]], self.dfcart[:,axes[2]], zdir='z', s=size, c='k', depthshade=True)
            ax.scatter(self.dfcart[ilabel,axes[0]], self.dfcart[ilabel,axes[1]], self.dfcart[ilabel,axes[2]], zdir='z', s=size*30, c='r', depthshade=True)
            ax.set_xlabel(self.data_name_cart[axes[0]], fontsize=35); ax.set_ylabel(self.data_name_cart[axes[1]], fontsize=25); ax.set_zlabel(self.data_name_cart[axes[2]], fontsize=25)
        
        ax.set_title(self.name, fontsize=30)
        #ax.view_init(30, 60)
        plt.show()


    ##############################################
    def dbscan_labels(self, eps=0.15, min_samples=15, display=True) :
        "Compte a DBSCAN clustering and return the largest cluster found"
        
        db = cluster.DBSCAN(eps=eps, min_samples=min_samples).fit(self.dfnorm)
        labels = db.labels_
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        max_size = 0
        if n_clusters_ > 0 :
            for i in range(n_clusters_) :
                ilabel = np.where(labels == i)[0]
                label_size = len(ilabel)
                if label_size > max_size :
                    ilabel_final = np.copy(ilabel)
                    max_size = label_size
            if display : print(str(n_clusters_)+" clusters, size of the largest: "+str(len(ilabel_final)))
            return ilabel_final
        else :
            print("ERROR 0 cluster found with eps="+str(eps)+" and min_samples="+str(min_samples))
            return []





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
    
    
    
        
        
        
    
