## Methods and types to help with the GAIA data.
##
## Julia 1.0 compliant
## 
## Use the astropy python package to read a votable.

module gaiaClustering

using PyCall
using DataFrames
using Clustering
using Statistics

include("data.jl")

export read_votable , filter_data , add_cartesian , normalization_PerBlock , toto1

end