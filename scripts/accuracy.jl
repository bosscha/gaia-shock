
## standalone script to extract recursively the stellar clusters in a votable with or w/o
## optimization of the DBSCAN parameters.
##

using DataFrames, Query
using CSV, Glob, Dates
using Statistics, Random, UUIDs
using Printf, ArgParse
using Distributions

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
        "--pca", "-p"
            help = "Add Principal Component coordinates in oc file and save PC vectors"
            action = :store_true
        "--zpt", "-z"
            help = "Apply Zero Point offset correction (Lindengren 2020)"
            action = :store_true
        "--iso", "-i"
                help = "Isochrone fitting"
                action = :store_true    
        "--tail", "-t"
                help = "Apply a second step to extract more isolated members (tails e.g.)"
                action = :store_true                   
        "--maxdist", "-d"
            help = "maximum distance for stars in pc"
            arg_type = Float64
        "--mindist"
            help = "miniimum distance for stars in pc"
            arg_type = Float64
        "--qmetric", "-q"
            help = "Q metric method to choose best solution in final DBSCAN clusters. The options are: Qc, Qn, QcQn, QcQnhigh"
            arg_type = String
        "--qc"
        help = "minimum Qc for the optimisation (default=2.7)"
        arg_type = Float64
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
            required = false
    end

    return parse_args(s)
end
##
## Main function with the errors on parallax
##


function simulate_parallax_noise(true_parallax::Float64, gaia_mag::Float64; 
    num_observations::Int=100,
    base_parallax_error::Float64=0.03,  # Base error in milliarcseconds (mas)
    mag_dependence_factor::Float64=0.1, # Controls how much magnitude impacts error
    scan_angle_effect::Bool=true,        # Include scan-angle-dependent error
    scan_angle_std::Float64=0.015,    # Standard deviation of scan angle effect (mas)
    systematic_error_floor::Float64 = 0.01, # Minimum error (mas), even for bright stars
    add_outliers::Bool=false,          # Simulate occasional large errors
    outlier_probability::Float64=0.01,  # Probability of an outlier
    outlier_std_factor::Float64=10.0,   # Outlier standard deviation multiplier
    rng::AbstractRNG = Random.GLOBAL_RNG #allow passing a random number generator
    )

"""
Simulates parallax measurements with noise, mimicking Gaia data characteristics.

Args:
true_parallax: The true parallax of the star in milliarcseconds (mas).
gaia_mag: The Gaia G-band magnitude of the star.
num_observations: The number of simulated parallax observations.
base_parallax_error: The base parallax error for a bright star (mas).
mag_dependence_factor:  A factor controlling the magnitude dependence of the error.  
   Higher values mean magnitude has a larger impact.
scan_angle_effect: Whether to include a scan-angle-dependent error component.
scan_angle_std: Standard deviation of the scan angle effect (mas).
systematic_error_floor: Minimum error achievable, even for very bright stars (mas).
add_outliers:  Whether to include occasional large outlier measurements.
outlier_probability: Probability of a single observation being an outlier.
outlier_std_factor:  How many times larger the outlier standard deviation is 
 compared to the normal error.
rng:  Random number generator.  Defaults to the global RNG.  Pass a specific
instance (e.g., `MersenneTwister(1234)`) for reproducible results.


Returns:
A vector of simulated parallax measurements (mas).  The length of the vector 
will be equal to `num_observations`.

Example:
```julia
simulated_parallaxes = simulate_parallax_noise(1.0, 15.0)  # True parallax = 1 mas, G=15 mag
println("Mean simulated parallax: ", mean(simulated_parallaxes))
println("Standard deviation: ", std(simulated_parallaxes))

# Reproducible example with a specific seed:
using Random
rng = MersenneTwister(42);  # Use a specific seed
simulated_parallaxes_reproducible = simulate_parallax_noise(1.0, 15.0; rng=rng)
```
"""

# debug_red("### simulated noise ...........")
# Magnitude-dependent error.  Brighter stars have smaller errors.
mag_error = base_parallax_error * (1.0 + mag_dependence_factor * (gaia_mag - 12.0))
mag_error = max(mag_error, systematic_error_floor) # Apply error floor


# Generate the simulated measurements
measurements = Vector{Float64}(undef, num_observations)
for i in 1:num_observations
# Basic parallax error (normally distributed)
parallax_error = rand(rng, Normal(0.0, mag_error))

# Scan angle effect (if enabled)
if scan_angle_effect
scan_angle_error = rand(rng, Normal(0.0, scan_angle_std))
parallax_error += scan_angle_error
end

# Add outliers
if add_outliers && rand(rng) < outlier_probability
parallax_error += rand(rng, Normal(0.0, mag_error * outlier_std_factor))
end

measurements[i] = true_parallax + parallax_error
end

return measurements
end



#########
function _extra(m::GaiaClustering.meta, optim)
    println("## Extra version to estimate the effects of noise on parallax...")
    tstart= now()
    rng = MersenneTwister()
    uuid=uuid4(rng)
    m.uuid= uuid

    println("###########################")
    println("## Starting with $(m.votname)")
    println("## Starting at $tstart")
    println("## Id $uuid")

    df , dfcart , dfcartnorm = _get_data(m)

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

function _get_data(m::GaiaClustering.meta)
    println("##  to add errors on parallax measurements...")
    println("## Distance cut : $(m.mindist) $(m.maxdist) pc")

    if m.zpt=="yes"
        zoff= true
        println("## Applying Zero Point offset correction on parallax...")
    else
        zoff= false
    end

    data       = read_votable(m.votdir*"/"*m.votname)
    ngaia = length(data)
    debug_red("# data length  $ngaia")

    ##
    datanoise= data
    println("## Adding noise on parallax...")
    # datanoise.parallax = datanoise.parallax .+ simulate_parallax_noise(datanoise.parallax, datanoise.gaia_mag)
    # parallax[i] = convert(Float64, get(gaia,i-1).parallax)
    for i in 1:ngaia
        parallax = convert(Float64, get(datanoise,i-1).parallax) 
        magG = convert(Float64, get(datanoise,i-1).phot_g_mean_mag)
        # debug_red("--- $parallax $magG")
        if isnan(magG)
            println("#### NaN")
            magG= 20.0
        end
        parallaxnoise= simulate_parallax_noise(parallax, magG, num_observations=1)
        # debug_red(parallaxnoise)
        
        # a0= datanoise[i].parallax
        # a1= datanoise[i-1].parallax
        # println("$parallax,$a0,$a1")

        datanoise[i].parallax = parallaxnoise
    end
    ##

    df         = filter_data(datanoise, [m.mindist, m.maxdist],zpt=zoff)
    dfcart     = add_cartesian(df)
    blck       = [[1,2,3],[4,5], [6,7,8]]
    wghtblck   = [4.0,5.0,1.0]
    norm       = "identity"

    dfcartnorm , scale8 = normalization_PerBlock(dfcart, blck, wghtblck , norm, false)

    return(df, dfcart , dfcartnorm)
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
    qc= parsed_args["qc"]

    if parsed_args["o"] opt="yes" else opt= "no" end
    if parsed_args["pca"] pca="yes" else pca= "no" end
    if parsed_args["zpt"] zpt="yes" else zpt= "no" end
    if parsed_args["iso"] iso="yes" else iso= "no" end
    if parsed_args["tail"] tail="yes" else tail= "no" end

    ############
    header_extract()

    if metafile != nothing
        m= read_params(metafile, false)

        opt= m.optim
        pca= m.pca
        zpt= m.zpt
        iso= m.iso
        tail=m.tail
        if votable == nothing votable= m.votname end
    else
        println("## The default options are used.")
        m= set_default_params()
    end

    m.optim= opt
    m.votname= votable
    m.pca= pca ; m.zpt= zpt ; m.iso= iso ; m.tail= tail

    if m.optim == "yes"
        isoptimized= true
    else
        isoptimized= false
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
    if qc != nothing m.minQc= qc end   


    ## Simulations ...
    NSimul= 1000  # number of simulations
    for i in 1:NSimul
        println("## Simulation : $i")
        _extra(m, isoptimized)
        debug_red("Extra version to simulate and compute noise...")
    end
end
