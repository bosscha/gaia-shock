"""
Script to run tests on GAIA data.
"""


import gaia_utils as gu
import os

## directory
rootdir = "/home/stephane/Science/GAIA"
wdir    = "%s/products"%(rootdir)
datadir = "%s/master/notebooks/data"%(rootdir)

os.chdir(wdir)


###################################
clustername = "NGC 2682"
voname = 'NGC 2682-1.0deg.vot'
RADIUS   = 1.0
kCluster = 6
votable_disk = False
BINSIZE = 64
SIGMA = 1.0
distclust = 830.0
XYRANGE = [-20., 20]
WEIGHT = [20.,5.,5.,6.,6.]

################################

## Read the data and do the conversion


source = gu.source(clustername)
source.query(RADIUS, errtol = 0.5, dump = True)
source.read_votable(voname)
source.convert_filter_data()
source.normalization0_1()
