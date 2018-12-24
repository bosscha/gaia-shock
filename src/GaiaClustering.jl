## Methods and types about stellar clusters.
## 
##
## Julia 1.0 compliant
## 
## Use the astropy python package to read a votable.

module GaiaClustering

using PyCall
using DataFrames , Formatting

using Clustering
using Statistics , Distributions ,Random

import PyPlot

## include all the types
include("types.jl")

## GAIA function to deal with data
include("data.jl")
export read_votable , filter_data , add_cartesian , normalization_PerBlock , copy1 ,
  subsetDf 

## Geometry functions (Voronoi, correlation2d)
include("geometry.jl")
export voronoi

## stelllar cluster analysis
include("stellarcluster.jl")
export metric , clusters , find_clusters , find_cluster_label, get_properties_SC , metric2

## MCMC for gaia
include("mcmc.jl")
export theta , thetaiter , abc_mcmc_dbscan, ministats , abc_mcmc_dbscan_full , ministats_full , 
theta_full , thetaiter_full , getDfcartnorm

## plotting functions
include("plots.jl")
export show_text , plot_dbscan_mcmc , plot_cluster , plot_dbscanfull_mcmc

## utils methods
include("utils.jl")
export isnotnan

end