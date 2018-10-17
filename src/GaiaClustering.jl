## Methods and types about stellar clusters.
## 
##
## Julia 1.0 compliant
## 
## Use the astropy python package to read a votable.

module GaiaClustering

using PyCall
using DataFrames
using Clustering
using Statistics
using Distributions
using Random

## include all the types
include("types.jl")

## GAIA function to deal with data
include("data.jl")
export read_votable , filter_data , add_cartesian , normalization_PerBlock , copy1

## Geometry functions (Voronoi, correlation2d)
include("geometry.jl")
export voronoi

## stelllar cluster analysis
include("stellarcluster.jl")
export metric , clusters , find_clusters

## MCMC for gaia
include("mcmc.jl")
export theta , thetaiter , abc_mcmc_dbscan, ministats

end