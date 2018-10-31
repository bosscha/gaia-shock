************
Source
************

The MWSC list of stellar clusters are taken from HEASARCS. The maximum distance was set to 2 kpc to get sufficient acceptable GAIA data. 
In the **sample folder** the following nb:
* getVotableGaia.ipynb : python nb to loop over the list to download the VO table in a conesearch of 5 times the cluster radius or at least 1 degree. It will resume if interrupted.
* upgrade-SClist.ipynb : julia nb to upgrade the mwsc list adding the column voname to join with the results from the analysis.

Requirements
============
Python:
astroquery

Julia:
CSV
