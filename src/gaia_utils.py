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
        
        minRA = min(self.data[ "ra"])
        maxRA = max(self.data[ "ra"])
        minDec = min(self.data[ "dec"])
        maxDec = max(self.data[ "dec"])
        
        centerRA  = self.data[ "ra"].mean()
        centerDec = self.data[ "dec"].mean()
        
        halfRA  = minRA + (maxRA-minRA)/2.
        halfDec = minDec + (maxDec-minDec)/2.
        
        if abs(halfRA - centerRA ) < tol*math.cos(centerDec*DEG2RAD) and abs(halfDec - centerDec ) < tol:
            return(True)
        else:
            return(False)
        
        
        
        
    