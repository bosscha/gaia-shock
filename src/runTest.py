"""
Script to run tests on GAIA data.
"""


import gaia_utils as gu
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import pandas as pd

from astroquery.gaia import Gaia

c = coord.SkyCoord.from_name("NGC 2232")
radius = 0.2

queryaql = "SELECT * FROM gaiadr2.gaia_source WHERE \
CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS' \
,%f ,%f,%f))=1;"%(c.ra.deg+radius/10.,c.dec.deg+radius/10.,radius)

job = Gaia.launch_job_async( queryaql, dump_to_file=False)
# print (job)

data = job.get_results()


gaia = gu.gaiaSet(data)
im = gaia.sampling_filtering(data['l'],data['b'], 256, 3.0)


print(gaia.isHomogeneous())
