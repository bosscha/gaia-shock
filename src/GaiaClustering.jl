## Methods and types about stellar clusters.
##
## Julia 1.0 compliant
##
## Use the astropy python package to read a votable.

module GaiaClustering

using PyCall
using DataFrames , Formatting, Printf

using Clustering
using Statistics , Distributions ,Random
using Distances, LightGraphs, SimpleWeightedGraphs

## For wavelet in imaging
import Interpolations

using Images , Base

import PyPlot, StatsBase , CSV
import Distances

## include all the types
include("types.jl")
export modelCauchy , Grav , abcfull , mcfull, modelfull

## GAIA function to deal with data
include("data.jl")
export read_votable , filter_data , add_cartesian , normalization_PerBlock , copy1 ,
  subsetDf , galXYZ , PM_equatorial2galactic , galUVW ,export_df

## Geometry functions (Voronoi, correlation2d)
include("geometry.jl")
export voronoi

## stelllar cluster analysis
include("stellarcluster.jl")
export metric , clusters , find_clusters , find_cluster_label, get_properties_SC , metric2

## MCMC for gaia
include("mcmc.jl")
export theta , thetaiter , abc_mcmc_dbscan, ministats , abc_mcmc_dbscan_full , ministats_full ,
theta_full , thetaiter_full , getDfcartnorm , check_qminqstar_full

## plotting functions
include("plots.jl")
export show_text , plot_dbscan_mcmc , plot_cluster , plot_dbscanfull_mcmc

## imaging functions
include("imaging.jl")
export atrous , addWav , thresholdingWav , noiseWav

## utils methods
include("utils.jl")
export isnotnan , read_blacklist

## Mass segregation and stellar clustering
include("massSegregation.jl")
export mst_graph, lambda_mst, sample_size, kappa_ms, select_massivestars, get_kappaMS ,
  get_Q

 ## Spatial structure parameters
include("spatialStructure.jl")
export spatialParameter, density2D , density3D , locdens2d , locdens3d, model_rad, fdens1

## graph methods
include("graph.jl")
export zagreb_first, zagreb_second

end
