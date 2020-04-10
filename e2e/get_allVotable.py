#!/usr/bin/env python
# coding: utf-8

## End-to-End scripts

# GAIA Clusters <br>
# Download the votables by taking the RADIUS = 5. * radius_clusterc <br>
# The relative error on parallax and pm are 20% <br>

# In[ ]:

import sys, os
rootdir = os.environ.get("GAIA_ROOT")

sys.path.append("%s/master/src"%(rootdir))

import astropy.coordinates as coord
from astropy.io.votable import parse
from astropy.table import Table
from astropy import units as u
from astroquery.gaia import Gaia

import numpy as np
import pandas as pd
import shutil

import gaia_utils as gu


## directory
wdir    = "%s/e2e_products/votable-TEST"%(rootdir)
datadir = "%s/e2e_products"%(rootdir)

os.chdir(wdir)

#### Input Parameters########################################################
#### from check_sampleOCdr2.py
# filelist = datadir + "/"+"BrowseTargets.18292.1530479692.gaia.selected.txt"
# fileoutGaia = datadir + "/"+"BrowseTargets.18292.1530479692.gaia.votable.txt"
filelist = datadir + "/"+"sc-list-2020.csv"
fileoutGaia = datadir + "/"+"sc-list-2020.gaia.votable.csv"

RADIUS   = 1.0

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

###############################################################
######  Main loop
###############################################################

df_cluster = read_cluster_list(filelist)
print(df_cluster.index)
print(df_cluster.columns)

lastSC , lastrow = find_lastSC(fileoutGaia, filelist, df_cluster)

print(lastSC)
totalSelected= 0
totalDiscarded= 0

for index, row in df_cluster.iloc[lastrow:].iterrows():
    clustername = row['name'].strip()
    print("\n\n")
    print("Cluster: %s"%(clustername))
    print("Distance: %3.1f pc"%(row['distance']))

    radius = max(10 * float(row['cluster_radius']) , 2.0)   ## minimum  2 degree
    radius = min(45.0, radius)                              ## maximum 45 degree
    print("Field radius: %3.1f deg"%(radius))
    rasplit = row['ra'].split(' ')
    decsplit = row['dec'].split(' ')
    racluster = "%sh%sm%ss"%(rasplit[0],rasplit[1],rasplit[2])
    deccluster = "%sd%sm"%(decsplit[0],decsplit[1])
    c = coord.SkyCoord(racluster, deccluster, frame='icrs')

    cluster = gu.source(clustername)

    filename = cluster.query(radius, coordCluster = [c.ra.deg,c.dec.deg], errtol = 0.2, dump = True)
    gaia = gu.gaiaSet(cluster.data)
    selected = gaia.isHomogeneous(tol = 0.1)

    if selected:
        totalSelected += 1
        filedst = "%s-%3.1fdeg.vot"%(clustername, radius)
        shutil.move(filename,filedst)
        write_SCgaia(fileoutGaia,row)
        print("## Selected")
    else:
        totalDiscarded += 1
        os.remove(filename)
        print("## Not selected")

    print("## Total Selected %d - Discarded %d "%(totalSelected, totalDiscarded))
