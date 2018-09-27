## Methods and types about stellar clusters.
## 
##
## Julia 1.0 compliant
## 
## Use the astropy python package to read a votable.

module gaiaClustering

using PyCall
using DataFrames
using Clustering
using Statistics

## GAIA function to deal with data
include("data.jl")

export read_votable , filter_data , add_cartesian , normalization_PerBlock , copy

end