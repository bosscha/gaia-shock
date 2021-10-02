## wrapper and line command script to get GAIA votable in a conesearch

import argparse
import os, sys
import astropy.coordinates as coord
from astropy import units as u

rootdir = os.getenv('GAIA_ROOT')
sys.path.append("%s/run/src" % (rootdir))

import gaia_basics as gb

def get_options():
    parser = argparse.ArgumentParser(
        description='Wrapper to conesearch GAIA data and return a votable')
    parser.add_argument("-f", help="fieldname", default="gaia", type=ascii)
    parser.add_argument("-l", type=float)
    parser.add_argument("-b", type=float)
    parser.add_argument("--ra", type=float)
    parser.add_argument("--dec", type=float)
    parser.add_argument("--fov", default= 2.0, type=float)
    parser.add_argument("--tol", default= 0.2, type=float)

    args = parser.parse_args()
    return(args)


def coord_galactic(longi,lati):
    c = coord.SkyCoord(frame="galactic", l=longi*u.degree, b=lati*u.degree)
    return(c.fk5)

########################################################
if __name__ == "__main__":

    args= get_options()
    print(args.f)

    c=coord_galactic(args.l,args.b)
    
    print("%s--"%(args.f))

    cluster = gb.source(args.f)
    filename = cluster.query(args.fov, coordCluster = [c.ra.deg,c.dec.deg], errtol = args.tol, dump = True)
