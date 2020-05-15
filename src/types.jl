## definition of the different types and constant


DEG2RAD =  π / 180.
Grav    = 4.30091e-3          ## Gravitational constant in km/s, pc, Msun

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
    pml::Float64
    pmb::Float64
    vl::Float64
    vb::Float64
    vrad::Float64
    xdisp::Float64
    ydisp::Float64
    zdisp::Float64
    vldisp::Float64
    vbdisp::Float64
    vraddisp::Float64
end

## 2MASS  type
mutable struct twomass
    J::Array{Float64}
    Jerr::Array{Float64}
    H::Array{Float64}
    Herr::Array{Float64}
    K::Array{Float64}
    Kerr::Array{Float64}
end

### Spatial Structure Types
##
## position from the gravity centered
struct sc2dcentered
    nxy::Int
    xx::Array{Float64}
    yy::Array{Float64}
    radius::Array{Float64}
    dens::Array{Float64}
    nbin::Int
    radbin::Array{Float64}
    densbin::Array{Float64}
    densbinerr::Array{Float64}
end

## prior model / mcmc
mutable struct mcmcCauchy
    Cmean::Float64
    Cdisp::Float64
    smean::Float64
    sdisp::Float64
    mmean::Float64
    mdisp::Float64
    nburnin::Int
    niter::Int
end

mutable struct modelCauchy
    C::Float64
    s::Float64
    m::Float64
end

mutable struct mcCauchy
    C::Array{Float64}
    s::Array{Float64}
    m::Array{Float64}
end


### METADATA for settings, optimization, etc
mutable struct meta
    rootdir::String
    wdir::String
    votdir::String
    plotdir::String
    ocdir::String
    votname::String

    minQc::Float64
    minQn::Int
    maxQn::Int
    forcedminstars::Int
    mingoodsolution::Int
    niterqminq::Int
    nburnout::Int
    nchain::Int
    maxiter::Int

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

    aperture2d::Float64
    maxaperture2d::Float64
    aperture3d::Float64
    maxaperture3d::Float64
    aperturev::Float64
    maxaperturev::Float64
    nboot::Int
    labels::String

    meta()= new()
end

## Basic properties of a SC v2
struct SCproperties2
    nstars::Int
    distance::Float64
    ra::Float64
    dec::Float64
    l::Float64
    b::Float64
    parallax::Float64
    pmra::Float64
    pmdec::Float64
    pml::Float64
    pmb::Float64
    vl::Float64
    vb::Float64
    vrad::Float64
    xdisp::Float64
    ydisp::Float64
    zdisp::Float64
    vldisp::Float64
    vbdisp::Float64
    vraddisp::Float64
    offdeg::Float64      ## offset in degree wrt. field center

    SCproperties2()= new()
end
