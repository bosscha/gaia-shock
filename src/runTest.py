"""
Script to run tests on GAIA data.
"""


import gaia_utils as gu
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import pandas as pd
from gaia.tap import cone_search


c = coord.SkyCoord.from_name("M67")
data = cone_search(c.ra.deg, c.dec.deg, 0.5, table="gaiadr2.gaia_source", dump_to_file=True)

gaia = gu.gaiaSet(data)

print(gaia.isHomogeneous())
