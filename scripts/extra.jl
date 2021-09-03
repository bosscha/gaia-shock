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


##
## parse extra options/flags
##
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
##
## get the data. The weightings are fake. Should be applied later
##
function getdata(filevot)
    voname = filevot

    data       = read_votable(voname)
    df         = filter_data(data, [0,2100])
    dfcart     = add_cartesian(df)
    blck       = [[1,2,3],[4,5], [6,7,8]]
    wghtblck   = [4.0,5.0,1.0]
    norm       = "identity"

    dfcartnorm , scale8 = normalization_PerBlock(dfcart, blck, wghtblck , norm, false)
    return(df, dfcart , dfcartnorm)
end
##
## Main function
##
function main(m::GaiaClustering.meta, optim)
    tstart= now()
    println("###########################")
    println("## Starting with $(m.votname)")
    println("## Starting at $tstart")

    df , dfcart , dfcartnorm = getdata(m.votname)
    cycle, flag= cycle_extraction_optim(df, dfcart, m, optim)

    tend= now()
    println("## Ending at $tend")
    println("## number of cycle: $cycle , flag:$flag ")
    println("##########################")
    println("##")

    duration= Dates.value(tend-tstart) / (1000*3600)
    durationstr= @sprintf("%3.3f", duration)
    @printf("## %s \n",specialstr("Duration: $durationstr hours","YELLOW"))
    @printf("## %s \n",specialstr("Votable done: $(m.votname)","YELLOW"))
    println("##\n##")
end
##################################### MAIN
let
    parsed_args = parse_commandline()

    metafile= parsed_args["m"]
    votable= parsed_args["votable"]
    ncycle= parsed_args["n"]
    if parsed_args["o"] opt="yes" else opt= "no" end

###############
    header_extract()

    m= read_params(metafile, false)

    m.optim= opt
    m.cyclemax= ncycle
    m.votname= votable

    main(m, parsed_args["o"])
end
