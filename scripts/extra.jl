## standalone script to extract recursively the stellar clusters in a votable with or w/o
## optimization of the DBSCAN parameters.
##


using DataFrames, Query
using CSV, Glob, Dates
using Statistics, Random
using Printf, ArgParse

rootdir =  ENV["GAIA_ROOT"]

push!(LOAD_PATH,"$rootdir/run/src")
using GaiaClustering



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-m"
            help = "configuration file"
            arg_type = String
            default = "configAll.ext"
        "-n"
            help = "maximum number of cycles"
            arg_type = Int
            default = 2
        "-o"
            help = "optimization of the weightings/DBSCAN"
            action = :store_true
        "votable"
        help = "votable"
        required = true
    end

    return parse_args(s)
end

##################################### MAIN
let
    parsed_args = parse_commandline()

    metafile= parsed_args["m"]
    votable= parsed_args["votable"]
    ncycle= parsed_args["n"]
    if parsed_args["o"] optim="yes" else optim= "no" end

###############
    header_extract()
    m= read_params(metafile, true)

    m.optim= optim
    m.cyclemax= ncycle
    m.votname= votable


    # main(votsublist, metafile, prefile)
end
