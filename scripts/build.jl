## standalone script to (re)build sets of candidate cluster (cycle) from different
## fields
##

rootdir =  ENV["GAIA_ROOT"]

push!(LOAD_PATH,"$rootdir/run/src")
using GaiaClustering

push!(LOAD_PATH,"$rootdir/run/scripts")
using extra

println("### testing build...")
extra(0,false)
