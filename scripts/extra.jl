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
        "--ncycle", "-n"
            help = "maximum number of cycles"
            arg_type = Int
        "-o"
            help = "optimization of the weightings/DBSCAN"
            action = :store_true
        "--maxdist" , "-d"
            help = "maximum distance for stars in pc"
            arg_type = Float64
        "--mindist"
            help = "miniimum distance for stars in pc"
            arg_type = Float64
        "--qmetric", "-q"
            help = "Q metric method to choose best solution in final DBSCAN clusters. The options are: Qc, Qn, QcQn, QcQnhigh"
            arg_type = String
        "--w3d"
            help = "XYZ weighting"
            arg_type = Float64
        "--wvel"
            help = "Velocity weighting"
            arg_type = Float64
        "--whrd"
            help = "Magnitud/color weighting"
            arg_type = Float64
        "--eps"
            help = "Epsilon parameter of DBSCAN"
            arg_type = Float64
        "--mcl"
            help = "Min_cluster parameter of DBSCAN"
            arg_type = Int
        "--mnei"
            help = "Min_neighbor parameter of DBSCAN"
            arg_type = Int
        "votable"
            help = "votable"
            required = true
    end

    return parse_args(s)
end
##
## get the data. The weightings are fake. Should be applied later
##
function getdata(m::GaiaClustering.meta)
    println("## Distance cut: $(m.mindist) $(m.maxdist) pc")

    data       = read_votable(m.votname)
    df         = filter_data(data, [m.mindist, m.maxdist])
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

    df , dfcart , dfcartnorm = getdata(m)

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
    ncycle= parsed_args["ncycle"]
    qmetric= parsed_args["qmetric"]
    maxdist= parsed_args["maxdist"]
    mindist= parsed_args["mindist"]
    eps= parsed_args["eps"]
    mcl= parsed_args["mcl"]
    mnei= parsed_args["mnei"]
    w3d= parsed_args["w3d"]
    wvel= parsed_args["wvel"]
    whrd= parsed_args["whrd"]



    if parsed_args["o"] opt="yes" else opt= "no" end

    ############
    header_extract()

    if metafile != nothing
        m= read_params(metafile, false)
        opt= m.optim
    else
        println("## The default options are used.")
        m= set_default_params()
    end

    m.optim= opt
    m.votname= votable

    if m.optim == "yes"
        isoptimize= true
    else
        isoptimize= false
    end

    if qmetric != nothing m.labels= qmetric end
    if ncycle != nothing m.cyclemax= ncycle end
    if maxdist != nothing m.maxdist= maxdist end
    if mindist != nothing m.mindist= mindist end
    if eps != nothing m.eps= eps end
    if mcl != nothing m.mcl= mcl end
    if mnei != nothing m.mnei= mnei end
    if w3d != nothing m.w3d= w3d end
    if wvel != nothing m.wvel= wvel end
    if whrd != nothing m.whrd= whrd end

    main(m, isoptimize)
end
