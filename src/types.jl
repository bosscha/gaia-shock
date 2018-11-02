## definition of the different 


DEG2RAD =  π / 180.

## Gaia data struct.
mutable struct Df
    ndata::Int32
    data::Array{Float64}
    raw::Array{Float64}
    err::Array{Float64}
end

## dbcan parameters
mutable struct model
    eps::Float64
    min_nei::Int
    min_cl::Int
end

## dbcan and weighting parameters
mutable struct modelfull
    eps::Float64
    min_nei::Int
    min_cl::Int
    w3d::Float64
    wvel::Float64
    whrd::Float64
end

#### mutable Markov Chains to get the stationary df for the ABC-MCMC/dbscan opt,
mutable struct mc
    eps::Array{Float64}
    mne::Array{Int32}
    mcl::Array{Int32}
    qc::Array{Float64}
    qn::Array{Int32}
end

#### mutable Markov Chains to get the stationary df for the ABC-MCMC/dbscan opt,
mutable struct mcfull
    eps::Array{Float64}
    mne::Array{Int32}
    mcl::Array{Int32}
    w3d::Array{Float64}
    wvel::Array{Float64}
    whrd::Array{Float64}
    qc::Array{Float64}
    qn::Array{Int32}
end

## parameters for the ABC MCMC optimization of the dbscan
struct abc
    minQ::Float64
    minstars::Int
    epsmean::Float64
    epsdisp::Float64
    min_nei::Int
    min_cl::Int
    ncoredisp::Int
    nburnout::Int
    niter::Int
end

## parameters for the ABC MCMC optimization of the dbscan
struct abcfull
    minQ::Float64
    minstars::Int
    forcedminstars::Int    
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
    nburnout::Int
    niter::Int
end


## Basic properties of a SC
struct SCproperties
    nstars::Int
    distance::Float64
    ra::Float64
    dec::Float64
    l::Float64
    b::Float64
    parallax::Float64
    pmra::Float64
    pmdec::Float64
    vra::Float64
    vdec::Float64
    vrad::Float64
    xdisp::Float64
    ydisp::Float64
    zdisp::Float64
    vradisp::Float64
    vdecdisp::Float64
    vraddisp::Float64
end