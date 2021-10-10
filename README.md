# gaia-shock

Pipeline to extract stellar clusters from GAIA data as votable. The method is using a DBSCAN-optimized clustering method in the 8-fold dimension with photometry and astrometry data. The optimization the clustering parameters is performed with an Approximate Bayesian Computation method (Weyant et al. 2013).
The extraction process is currently with a stable version using the Julia language for speed constraint.
Be aware that some stellar clusters are not physically real and a final cleaning can be done with some Jupyter NB but are still experimental. A visual inspection of the final stellar cluster is still recommended.


Documentation and information can be found in the wiki pages.
