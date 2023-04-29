## Methods and types about stellar clusters.
##
## Julia 1.0 compliant
##
## Use the astropy python package to read a votable.

module GaiaClustering

using DataFrames , Formatting, Printf , Dates ,  DelimitedFiles

using Clustering
using Statistics , Distributions ,Random, UUIDs
using Distances, LightGraphs, SimpleWeightedGraphs, NearestNeighbors
using Glob , GZip , FileIO , CSVFiles , ProgressMeter

## For wavelet in imaging
import Interpolations
using  Images, Base

using  StatsBase , CSV,  Distances , MultivariateStats , FHist

using PyCall , PyPlot

using VoronoiCells
import GeometryBasics as gb

VERSION= "1.7.5" 

## include all the types
include("types.jl")
export modelCauchy , Grav , abcfull , mcfull, modelfull, meta

## GAIA function to deal with data
include("data.jl")
export read_votable , filter_data , add_cartesian , normalization_PerBlock , copy1 ,
  subsetDf , galXYZ , PM_equatorial2galactic , galUVW ,export_df , equatorial2galactic ,
  angle4sphere , get_data

## Geometry functions (Voronoi, correlation2d)
include("geometry.jl")
export voronoi

## stelllar cluster analysis
include("stellarcluster.jl")
export metric , clusters , find_clusters2 , find_cluster_label, get_properties_SC , metric2 ,
find_cluster_label2, get_metrics,  get_properties_SC2,  score_cycle ,
remove_stars , edge_ratio, save_cycle , compute_PC , cycle_extraction_optim, save_cycle_optim , distance_cluster

## MCMC for gaia
include("mcmc.jl")
export theta , thetaiter , abc_mcmc_dbscan, ministats  , ministats_full ,
theta_full , thetaiter_full , getDfcartnorm , abc_mcmc_dbscan_full2,
check_qminqstar_full2, create_DFchain, extraction_mcmc

## plotting functions
include("plots.jl")
export show_text , plot_dbscan_mcmc , plot_cluster , plot_dbscanfull_mcmc , plot_cluster2 ,
plot_rawdata, plot_astrom, plot_tail , level_dens , plot_sky

## imaging functions
include("imaging.jl")
export atrous , addWav , thresholdingWav , noiseWav

## utils methods
include("utils.jl")
export isnotnan , read_blacklist , convertStruct2Df , specialstr , bold , yellow, purple, cyan, red, blue, header_extract, debug_red, checkdir

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

## tail methods
include("tail.jl")
export tail_stars , transform_df , distance_cmd_tail , surface_density

## functions for isochrone fitting
include("isochrones.jl")
export mist_df ,  read_isochrones, update_mag , weight_cmd, dist_cmd2iso , fit_isochrone , read_serial_mist , update_nan_oc , get_star_mass ,  perform_isochrone_fitting

## functions to be used in the build script mainly
include("_build.jl")
export extra , get_gaia_data , get_gaia_data_many , get_random_field , galactic2equatorial , rm_duplicated, get_chunks

## deprecated functions
include("deprecated.jl")
export abc_mcmc_dbscan_full , check_qminqstar_full , find_clusters , cycle_extraction

## testing functions
include("testing.jl")
export __plot_check , __plot_nstars , __tail_stars , __density_count , __level_dens, __plot_surface_density , __plot_mcmc_article
end
