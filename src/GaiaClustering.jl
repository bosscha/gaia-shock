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
using Random

## GAIA function to deal with data
include("data.jl")
export read_votable , filter_data , add_cartesian , normalization_PerBlock , copy1

## Geometry functions (Voronoi, correlation2d)
include("geometry.jl")
export voronoi

## stelllar cluster analysis
include("stellarcluster.jl")
export metric , clusters


end