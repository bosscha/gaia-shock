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

#### mutable Markov Chains to get the stationary df for the ABC-MCMC/dbscan opt,
mutable struct mc
    eps::Array{Float64}
    mne::Array{Int32}
    mcl::Array{Int32}
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