"""

Methods to help with the GAIA data.
Python3 compl

HISTORY:
    14.06.2018:
        - first shot
        - adding a class for a dataset from gaia-on-tap
        -
        
"""

__author__  = "SL: ALMA"
__version__ = "0.0.1@2018.06.14"

import math

DEG2RAD = math.pi / 180.

class gaiaSet:
    
    def __init__(self, data):
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
        
        
        
        
    
