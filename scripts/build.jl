## standalone script to (re)build sets of candidate cluster (cycle) from different
## fields
##

using DataFrames , TOML , ArgParse

rootdir =  ENV["GAIA_ROOT"]

push!(LOAD_PATH,"$rootdir/run/src")
using GaiaClustering

#################################### MAIN 
let
    println(ARGS)
    println("####################")
    println("### testing build...")

    metabuild = TOML.parsefile(ARGS[1])

    println(metabuild)

    extra(0,false)
end
