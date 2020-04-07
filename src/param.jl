## functions to set/read the different parameters for the extraction
##


mutable struct extractParam
    rootdir::String
    wdir::String
    votdir::String
    plotdir::String
    ocdir::String

    minQ::Float64
    minstars::Int
    maxstars::Int
    forcedminstars::Int
    mingoodsolution::Int
    niterqminq::Int
    nburnout::Int
    niter::Int

    epsmean::Float64
    epsdisp::Float64
    min_nei::Int
    min_cl::Int
    ncoredisp::Int
    w3dmean::Float64
    w3ddisp::Float64
    wvelmean::Float64
    wveldisp::Float64
    whrdmean::Float64
    whrddisp::Float64
end

## return the default parameters
##
function set_default_params()
    ## Directories
    rootdir =  ENV["GAIA_ROOT"]
    wdir    = "$rootdir/e2e_products"
    votdir  = "$wdir/votable"
    plotdir = "$wdir/plotsSelect"
    ocdir   = "$wdir/ocfull"

    ## MCMC
    ##
    minQ    = 2.7
    minstars = 40
    maxstars = 10000
    forcedminstars = 30
    mingoodsolution = 10        # minimum good solution for Q in check_qminqstar
    niterqminq = 500            # check_qminqstar..

    nburnout  = 2000            # burn-in iterations
    niter     = 15000           # MCMC iterations

    ## prior settings
    ##
    epsmean   = 2.5             # epsilon mean
    epsdisp   = 1.5             # epsilon dispersion
    min_nei   = 10              # min_nei mean
    min_cl    = 15              # min_cl mean
    ncoredisp = 10              # dispersion for min_nei and min_cl
    w3dmean   = 6.0
    w3ddisp   = 4.0
    wvelmean  = 6.0
    wveldisp  = 4.0
    whrdmean  = 2.0
    whrddisp  = 1.5

    ## Metrics find_clusters
    aperture2d    = 1.5
    maxaperture2d = 15
    aperture3d    = 3
    maxaperture3d = 20
    aperturev     = 3.0
    maxaperturev  = 20
    nboot         = 50          # bootstrap number

end

###
## Read a .ext file with some of the parameters set (others are default)
function read_params(filename)
end
