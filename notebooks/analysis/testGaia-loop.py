
# coding: utf-8

# ### Tests

# In[1]:


import sys, os
sys.path.append('../../src')

import astropy.coordinates as coord

from astropy.io.votable import parse
from astropy.table import Table
from astropy import units as u

import matplotlib.pyplot as plt
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd

from math import ceil
import gaia_utils as gu
from sklearn.cluster import KMeans

# get_ipython().run_line_magic('matplotlib', 'inline')

## directory
rootdir = "/home/stephane/Science/GAIA"
wdir    = "%s/products"%(rootdir)
datadir = "%s/master/notebooks/data"%(rootdir)

os.chdir(wdir)

#### Input Parameters########################################################
filelist = datadir+"/"+"BrowseTargets.18292.1530479692.gaia.selected.txt"
fileoutGaia = wdir + "/"+"BrowseTargets.18292.1530479692.gaia.skmeansPlots.txt"

RADIUS   = 1.0
kCluster = 12
votable_disk = True


# In[2]:


from astroquery.gaia import Gaia

# tables = Gaia.load_tables(only_names=True)

#for table in (tables):
#    print (table.get_qualified_name())
    


# In[3]:


def get_data(ra, dec, radius):
    "Get the GAIA data, we offset of 10% the center.."
    
    queryaql = "SELECT * FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS', gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS' ,%f ,%f,%f))=1;"%(ra+radius/10.,dec+radius/10.,radius)  

    job = Gaia.launch_job_async( queryaql, dump_to_file=False)
    r = job.get_results()
    
    return(r)


# In[4]:


## read the cluster list from HEASARC
def read_cluster_list(filelist):
    
    df = pd.read_csv(filelist, sep='|')
    
    return(df)

######
def init_SCgaia(filegaia, filelist):
    "init output if not there.."
    
    if not os.path.exists(filegaia):
        with open(filelist,"r") as f:
            header = f.readline()
        with open(filegaia,"w") as f:
            f.write(header)
    
#######
## write in fileoutGaia the selected cluster for GAIA
def write_SCgaia(filegaia, row):
    df = pd.DataFrame(row).T
    df.to_csv(filegaia,sep ="|", mode = "a", header= False, index = False)
    
    
######       
## find the last SC found.
def find_lastSC(filegaia,filelist, dfsc):
    if os.path.exists(filegaia):
        print(filegaia)
        dfcurrent = read_cluster_list(filegaia)
        last_cluster = dfcurrent['name'].iloc[[-1]].iloc[0]
        index = dfsc.index[dfsc['name'].str.contains(last_cluster)]
    else:
        init_SCgaia(filegaia, filelist)
        last_cluster = "No cluster"
        index = [0]
        
    return(last_cluster, index[0])


# In[5]:


## Filter the data
##

def plot_check(a1, a2, ifinal):
    "check the filltering"
    
    
    rcParams['figure.figsize'] = 14, 14
    f, axarr = plt.subplots(2, 2)

# ax.scatter(cluster_candidates["ra"], cluster_candidates["dec"], s=1, c="#000000")
    print(len(a1))
    print(len(a2))
    
    axarr[0,0].scatter(a1, a2, s=1, c="#000000")
    axarr[0,0].set_xlabel("a1")
    axarr[0,0].set_ylabel("a2")
    axarr[0,0].set_ylim([0,2000])

    axarr[0,1].scatter(a1[ifinal], a2[ifinal], s=1, c="#000000")
    axarr[0,1].set_xlabel("a1[ifinal]")
    axarr[0,1].set_ylabel("a2[ifinal]")

    plt.show()
    
def filter_data(lgal, bgal, distance, vra, vdec, cartesian = False, dist_range = [0., 2000], vra_range = [-200,200], vdec_range = [-200.,200]):
    "filter the data applying the range and return the sklearn-centric array"
    
    i1 = np.where((distance >= dist_range[0]) & (distance < dist_range[1]))
    i2 = np.where((vra >= vra_range[0]) & (vra < vra_range[1]))
    i12 = np.intersect1d(i1,i2)
    i3 = np.where((vdec >= vdec_range[0]) & (vdec < vdec_range[1]))
    ifinal = np.intersect1d(i12,i3)
    
    # plot_check(lgal,vra,ifinal)
    
    datask = np.zeros((len(ifinal),5))
    
    if cartesian:
        xx, yy, zz = convert_to_cartesian(lgal[ifinal], bgal[ifinal], distance[ifinal])
        datask[:,0] = xx
        datask[:,1] = yy
        datask[:,2] = zz
        datask[:,3] = vra[ifinal]
        datask[:,4] = vdec[ifinal]
    else:
        datask[:,0] = lgal[ifinal]
        datask[:,1] = bgal[ifinal]
        datask[:,2] = distance[ifinal]
        datask[:,3] = vra[ifinal]
        datask[:,4] = vdec[ifinal]
    
    return(datask)
    
    


# In[6]:


## plot2D and plot3D


def plot2d_labels(a1,a2,labels,nclusterk,centroid,xylab = ["a1","b1"],title = "Title", clustername= "NGC"):
    figname = "%s-2Dkmeans-gaia.png"%(clustername)
    rcParams['figure.figsize'] = 14, 21
    nrow = int(ceil(nclust / 3))
    ncol = 3           
    f, axarr = plt.subplots(nrow, ncol)

    for i in range(nclusterk):
        ilabel = np.where(labels == i)[0]
        row = int(ceil(i / 3)) - 1
        col = i % 3
        axarr[row,col].scatter(a1[ilabel],a2[ilabel], s=1, c="#000000")
        axarr[row,col].set_xlabel(xylab[0])
        axarr[row,col].set_ylabel(xylab[1])
        
        txt = "Stars:%d Dist:%3.1f Vra:%3.1f Vdec:%3.1f"%(len(ilabel), centroid[i][2],centroid[i][3],centroid[i][4])
        axarr[row,col].text(-0.1,1.02, txt, size=12, ha="left", 
         transform=axarr[row,col].transAxes)
        
    plt.savefig(figname)
    plt.show()

    
def plot2d_filtered_labels(a1,a2,labels,nclusterk,centroid, xylab = ["a1","b1"],title = "Title", clustername= "NGC"):
    
    "filtering with wavelet."
    
    gaia = gu.gaiaSet()
    
    figname = "%s-2Dkmeans-filtered-gaia.png"%(clustername)
    rcParams['figure.figsize'] = 14, 21
    nrow = int(ceil(nclust / 3))
    ncol = 3           
    f, axarr = plt.subplots(nrow, ncol)

    for i in range(nclusterk):
        ilabel = np.where(labels == i)[0]
        row = int(ceil(i / 3)) - 1
        col = i % 3
        im , imf = gaia.sampling_filtering(a1[ilabel],a2[ilabel], 256, 3.0)
        
        axarr[row,col].imshow(imf, cmap = 'gist_stern')
        axarr[row,col].set_xlabel(xylab[0])
        axarr[row,col].set_ylabel(xylab[1])
        
        txt = "Stars:%d Dist:%3.1f Vra:%3.1f Vdec:%3.1f"%(len(ilabel), centroid[i][2],centroid[i][3],centroid[i][4])
        axarr[row,col].text(-0.1,1.02, txt, size=12, ha="left", 
         transform=axarr[row,col].transAxes)
        
    plt.savefig(figname)
    plt.show()
  
def plot3d_labels(a1,a2,a3,labels,nclusterk,centroid,xylab = ["X","Y","Z"],title = "Title", clustername= "NGC"):
    figname = "%s-3Dkmeans-gaia.png"%(clustername)
    rcParams['figure.figsize'] = 14, 21
    nrow = int(ceil(nclust / 3))
    ncol = 3
    
    fig = plt.figure()

    for i in range(nclusterk):
        ilabel = np.where(labels == i)[0]
        row = int(ceil(i / 3)) 
        col = i % 3 + 1
        axarr = fig.add_subplot(nrow,ncol,i+1, projection='3d')
        axarr.scatter(a1[ilabel],a2[ilabel],a3[ilabel], c="r",marker ="*")
        axarr.set_xlabel(xylab[0])
        axarr.set_ylabel(xylab[1])
        axarr.set_zlabel(xylab[2])
    
    plt.savefig(figname)
    plt.show()


# In[7]:


## astrometric conversion
## 
def convert_to_cartesian(lgal, bga, dist):
    "Convert ra,dec (ICRS) and distance (pc) to Cartesian reference"
    
    xx = np.zeros(len(lgal))
    yy = np.zeros(len(lgal))
    zz = np.zeros(len(lgal))
    
    for i in range(len(lgal)):
        c = coord.SkyCoord(l=lgal[i]*u.degree, b=bgal[i]*u.degree, distance=dist[i]*u.pc, frame='galactic')
        
        xx[i] = c.cartesian.x.value
        yy[i] = c.cartesian.y.value
        zz[i] = c.cartesian.z.value
        
    return(xx,yy,zz)


# In[ ]:


###############################################################
######  Main loop
###############################################################

df_cluster = read_cluster_list(filelist)
print(df_cluster.index)
print(df_cluster.columns)

lastSC , lastrow = find_lastSC(fileoutGaia, filelist, df_cluster)

print(lastSC)

for index, row in df_cluster.iloc[lastrow:].iterrows():
    clustername = row['name'].strip()
    print("Cluster: %s"%(clustername))
    print("Distance: %3.1f pc"%(row['distance']))
    
    rasplit = row['ra'].split(' ')
    decsplit = row['dec'].split(' ')
    racluster = "%sh%sm%ss"%(rasplit[0],rasplit[1],rasplit[2])
    deccluster = "%sd%sm"%(decsplit[0],decsplit[1])
    c = coord.SkyCoord(racluster, deccluster, frame='icrs')
    
    if not votable_disk:
        data = get_data(c.ra.deg, c.dec.deg,radius= RADIUS)
    else:
        voname = "%s-%3.1fdeg.vot"%(clustername,RADIUS)
        votable = parse(voname)
        for table in votable.iter_tables():
            data = table.array
                    
    lgal = data['l']
    bgal = data['b']
    pmas = data['parallax']
    distance = 1000. / np.ma.filled(pmas, -9999.)
    pmra = np.ma.filled(data['pmra'], -9999999.)
    pmdec= np.ma.filled(data['pmdec'],-9999999.)
    vdec = 4.74 * pmdec / pmas   ##?
    vra  = 4.74 * pmra  / pmas


    print("## Total stars: %d"%(len(lgal)))

    # datask = np.zeros((len(vra),5))
    datask = filter_data(lgal,bgal,distance,vra,vdec ,cartesian = False)
    print("## Stars selected: %d"%(len(datask[:,0])))


    ## fitted cluster in k-means
    nclust = kCluster

    print("## computing k-means...")
    kmeans = KMeans(n_clusters=nclust, max_iter = 500, n_init = 50)
    kmeans.fit(datask)

    centroid = kmeans.cluster_centers_
    labels = kmeans.labels_

    print("## Centroids:")
    print(centroid)
    # print(labels)
    for i in range(nclust):
        nstar = len(np.where(labels == i)[0])
        print("## Label: %d, %d stars"%(i,nstar))

    ## plot the different clusters

    a1 = datask[:,0]     ## l
    a2 = datask[:,1]     ## b
    plot2d_labels(a1,a2,labels,nclust,centroid,xylab = ["l","b"],title = "Title", clustername= "%s-%3.0fpc"%(clustername,row['distance']))
    plot2d_filtered_labels(a1,a2,labels,nclust,centroid,xylab = ["l","b"],title = "Title", clustername= "%s-%3.0fpc"%(clustername,row['distance']))

    write_SCgaia(fileoutGaia,row)
    

