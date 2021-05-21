#!/usr/bin/env python
# coding: utf-8

# ## Sample  GAIA DR2
#
# To search for specific data in DR2. <br>
# Note that gaia-on-tap is by default for DR1 <br>

import sys, os
rootdir = os.environ.get("GAIA_ROOT")

sys.path.append("%s/run/src"%(rootdir))

import astropy.coordinates as coord
from astropy.coordinates import Angle

import matplotlib.pyplot as plt
from pylab import rcParams
import pandas as pd
from gaia.tap import cone_search

import gaia_utils as gu

## directory
wdir    = "%s/e2e_products"%(rootdir)
datadir = "%s/master/notebooks/data"%(rootdir)

os.chdir(wdir)

# ### OC selection


def plotCluster(cluster_candidates, clustername, display = True, dirplot = ""):

    figname = "%s-gaia.png"%(clustername)

    rcParams['figure.figsize'] = 14, 21
    f, axarr = plt.subplots(3, 2)

# ax.scatter(cluster_candidates["ra"], cluster_candidates["dec"], s=1, c="#000000")
    axarr[0,0].scatter(cluster_candidates["ra"], cluster_candidates["dec"], s=1, c="#000000")
    axarr[0,0].set_xlabel(r"$\alpha$")
    axarr[0,0].set_ylabel(r"$\delta$")

    axarr[0,1].scatter(cluster_candidates["l"], cluster_candidates["b"], s=1, c="#000000")
    axarr[0,1].set_xlabel(r"l")
    axarr[0,1].set_ylabel(r"b")

    axarr[1,0].scatter(cluster_candidates["l"], cluster_candidates["parallax"], s=1, c="#000000")
    axarr[1,0].set_xlabel(r"l")
    axarr[1,0].set_ylabel(r"p (mas)")
    axarr[1,0].set_ylim([-1,4])

    axarr[1,1].scatter(cluster_candidates["ra"], cluster_candidates["pmra"], s=1, c="#000000")
    axarr[1,1].set_xlabel(r"$\alpha$")
    axarr[1,1].set_ylabel(r"PM RA (mas/yr)")
    axarr[1,1].set_ylim([-40,40])

    axarr[2,0].scatter(cluster_candidates["pmdec"], cluster_candidates["pmra"], s=1, c="#000000")
    axarr[2,0].set_xlabel(r"PM DEC (mas/yr)")
    axarr[2,0].set_ylabel(r"PM RA (mas/yr)")
    axarr[2,0].set_xlim([-40,40])
    axarr[2,0].set_ylim([-40,40])

    axarr[2,1].scatter(cluster_candidates["parallax"], cluster_candidates["pmra"], s=1, c="#000000")
    axarr[2,1].set_xlabel(r"p (mas)")
    axarr[2,1].set_ylabel(r"PM RA (mas/yr)")
    axarr[2,1].set_xlim([-1,4])
    axarr[2,1].set_ylim([-40,40])

    f.subplots_adjust(hspace=0.5)

    plt.savefig(wdir+dirplot+"/"+figname)
    if display:
        plt.show()

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
        dfcurrent = read_cluster_list(filegaia)
        last_cluster = dfcurrent['name'].iloc[[-1]].iloc[0]
        index = dfsc.index[dfsc['name'].str.contains(last_cluster)]
    else:
        init_SCgaia(filegaia, filelist)
        last_cluster = "No cluster"
        index = [0]

    return(last_cluster, index[0])

#####################################################
# ### Main ..
#
# Making the plots and analysis of the targets.

filelist = datadir+"/"+"BrowseTargets.18292.1530479692"

fileoutGaia = wdir + "/"+"BrowseTargets.18292.1530479692.gaia.selected.txt"

df_cluster = read_cluster_list(filelist)
print(df_cluster.index)
print(df_cluster.columns)

lastSC , lastrow = find_lastSC(fileoutGaia, filelist, df_cluster)

print(lastrow)

for index, row in df_cluster.iloc[lastrow:].iterrows():
    clustername = row['name'].strip()
    print("Cluster: %s"%(clustername))

    # c = coord.SkyCoord.from_name(clustername)
    rasplit = row['ra'].split(' ')
    decsplit = row['dec'].split(' ')
    racluster = "%sh%sm%ss"%(rasplit[0],rasplit[1],rasplit[2])
    deccluster = "%sd%sm"%(decsplit[0],decsplit[1])
    a = Angle(racluster)
    d = Angle(deccluster)
    print(a.degree)
    print(d.degree)

    sc  = gu.source(clustername)
    res = sc.query(0.5, coordCluster = [a.degree,d.degree],errtol = 0.2 )

    gaia = gu.gaiaSet(sc.data)
    selected = gaia.isHomogeneous(tol = 0.05)

    if selected:
        write_SCgaia(fileoutGaia,row)
        print("## Selected")
    else:
        print("## Not selected")

    plotCluster(sc.data, clustername, display = False ,dirplot="/plotsSelect")


    if len(sc.data) == 0:
        print("## No data...")

print("## GAIA selectopn done.")
