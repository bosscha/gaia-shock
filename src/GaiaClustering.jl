## Methods and types about stellar clusters.
##
## Julia 1.0 compliant
##
## Use the astropy python package to read a votable.

module GaiaClustering

using DataFrames , Formatting, Printf , Dates ,  DelimitedFiles

using Clustering
using Statistics , Distributions ,Random
using Distances, LightGraphs, SimpleWeightedGraphs

## For wavelet in imaging
import Interpolations
using  Images, Base

using  StatsBase , CSV,  Distances , MultivariateStats

using PyCall , PyPlot

VERSION= "1.4.3 "

## include all the types
include("types.jl")
export modelCauchy , Grav , abcfull , mcfull, modelfull, meta

## GAIA function to deal with data
include("data.jl")
export read_votable , filter_data , add_cartesian , normalization_PerBlock , copy1 ,
  subsetDf , galXYZ , PM_equatorial2galactic , galUVW ,export_df , equatorial2galactic ,
  angle4sphere

## Geometry functions (Voronoi, correlation2d)
include("geometry.jl")
export voronoi

## stelllar cluster analysis
include("stellarcluster.jl")
export metric , clusters , find_clusters2 , find_cluster_label, get_properties_SC , metric2 ,
find_cluster_label2, get_metrics,  get_properties_SC2, cycle_extraction, score_cycle ,
remove_stars , edge_ratio, save_cycle , compute_PC , cycle_extraction_optim, save_cycle_optim

## MCMC for gaia
include("mcmc.jl")
export theta , thetaiter , abc_mcmc_dbscan, ministats  , ministats_full ,
theta_full , thetaiter_full , getDfcartnorm , abc_mcmc_dbscan_full2,
check_qminqstar_full2, create_DFchain, extraction_mcmc

## plotting functions
include("plots.jl")
export show_text , plot_dbscan_mcmc , plot_cluster , plot_dbscanfull_mcmc , plot_cluster2 , plot_rawdata
## imaging functions
include("imaging.jl")
export atrous , addWav , thresholdingWav , noiseWav

## utils methods
include("utils.jl")
export isnotnan , read_blacklist , convertStruct2Df , specialstr , bold , yellow, purple, cyan, red,  header_extract

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

## metadata methods
include("metadata.jl")
export set_default_params , read_params , set_param!

## deprecated functions
include("deprecated.jl")
export abc_mcmc_dbscan_full , check_qminqstar_full , find_clusters

## testing functions
include("testing.jl")
export __plot_check , __plot_nstars
end
