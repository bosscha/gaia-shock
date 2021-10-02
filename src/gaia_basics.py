## basic module for getgaia extracter from gaia_utils

import astropy.coordinates as coord

from astropy.io.votable import parse
from astropy.table import Table
from astropy import units as u
from astroquery.gaia import Gaia
import numpy as np
import shutil


class source:
    "Class to download gaia data for a source and to perform a clustering extractopn"

    ########################
    def __init__(self, name):
        self.name = name
        ################################

    def query(self, radius, coordCluster=[], errtol=0.2, dump=False, table="gaiaedr3.gaia_source"):
        "do a conesearch"

        self.weighted = False

        if len(coordCluster) == 2:
            c = coord.SkyCoord(
                ra=coordCluster[0]*u.degree, dec=coordCluster[1]*u.degree, frame='icrs')
            c_icrs = c
        else:
            try:
                c = coord.SkyCoord.from_name(self.name)
                c_icrs = coord.SkyCoord(
                    ra=c.ra.deg * u.degree, dec=c.dec.deg * u.degree, frame='icrs')
            except:
                print("## Cluster name not found ...")

        pos = c_icrs.galactic
        self.l_cluster = pos.l.value
        self.b_cluster = pos.b.value
        self.ra = c.ra.deg
        self.dec = c.dec.deg

        queryaql = "SELECT * FROM {table} WHERE CONTAINS(POINT('ICRS',{table}.ra,{table}.dec),  \
                                    CIRCLE('ICRS',{ra:.10f},{dec:.10f},{radius:.10f})) = 1  AND abs(pmra_error/pmra)<{err:10f}  AND abs(pmdec_error/pmdec)< {err:.10f} AND abs(parallax_error/parallax)< {err:.10f};".format(table=table, ra=c.ra.deg, dec=c.dec.deg, radius=radius, err=errtol)

        print(queryaql)
        job = Gaia.launch_job_async(queryaql, dump_to_file=dump)
        self.data = job.get_results()

        if dump:
            filename = job.outputFile
            filedst = "%s-%3.1fdeg.vot" % (self.name, radius)
            shutil.move(filename, filedst)
            print("## %s created" % (filedst))
        else:
            filedst = "dummy.vot"
            filename = None

        print("## Query for %s done" % (self.name))
        print("## Total stars: %d" % (len(self.data)))

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

        print("## %s read..." % (voname))
        print("## Total stars: %d" % (len(data)))

        return(len(data))

    #####################################
