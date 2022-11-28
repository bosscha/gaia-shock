## functions to fit the 2/3D spatial structure of OCs
##

## extract parameters
##

function spatialParameter(ocfile ; ntest=10 , nbin=20, niter=10000 ,verbose=true, nburnin= 5000 , dfoc=-1, centered=false)
    let
        solC= []; sols= []; solm= []
        ηBest= 1e9 ; θBest= 0
        fitfound= false

        if dfoc == -1
            oc= CSV.read(ocfile, DataFrame ,  delim= ";")
            if verbose println("## file $ocfile read...") end
        else
            oc=dfoc
        end

        ## binned density
        r2d,ρ2d,err2d= density2D(oc.Y, oc.Z, nbin)    ## binned radius, default in log
        # r3d,ρ3d,err3d= density3D(oc.X, oc.Y, oc.Z, nbin)

        ## local density
        ## not used yet...
        nloc= 3
        rad_2d, loc_2d= locdens2d(oc.Y, oc.Z, nloc)    ## density estimated for each point
        # rad_3d, loc_3d= locdens3d(oc.X,oc.Y, oc.Z, nloc)

        oc2d= sc2dcentered(length(oc.Y), oc.Y, oc.Z, rad_2d, loc_2d , nbin , r2d, ρ2d, err2d )
   
        if verbose println("## Density profiles computed...") end

        ## prior estimation
        ##
        prior, priorDisp= priorGuess(r2d, ρ2d)

        if verbose
            println("## Prior Guess:")
            println(prior)
            println("## prior estimated...")
        end

        C_mean= prior.C ; C_disp= priorDisp.C
        s_mean= prior.s ; s_disp= priorDisp.s
        m_mean= prior.m ; m_disp= priorDisp.m

        pinit= mcmcCauchy(C_mean, C_disp, s_mean, s_disp, m_mean, m_disp, nburnin, niter)

        ## if true it takes the mean as first step, otherwise random??

        if verbose println("## Running mcmc 2D fitting") end

        for i in 1:ntest
            mci ,  θbest , ηbest = main_mcmc(oc2d, pinit, likelihood2dbin, false, verbose)

            lags= collect(1:niter)
            acfC= StatsBase.autocor(mci.C, lags)[end]
            acfs= StatsBase.autocor(mci.s, lags)[end]
            acfm= StatsBase.autocor(mci.m, lags)[end]

            if verbose
                println("## test: $i")
                println(θbest)
                println("## Likelihood: $ηbest")
                println("## ACF(C,s,m): $acfC , $acfs , $acfm")
            end

            ## AFC < 0.1 to be considered as not correlated...
            if abs(acfC) < 0.1 && abs(acfs) < 0.1 && abs(acfm) < 0.1 && θbest != 0
              push!(solC,abs(θbest.C))
              push!(sols,abs(θbest.s))
              push!(solm,abs(θbest.m))

                if ηbest < ηBest
                    fitfound= true
                    ηBest= ηbest
                    θBest= θbest
                end
            end
        end

        if !fitfound
          θBest= modelCauchy(1,1,1)
        end
        if length(solC) > 1
          θerr= modelCauchy(std(solC),std(sols),std(solm))
        else
          θerr= modelCauchy(1,1,1)
        end

        if θBest.C<0 || θBest.s<0 || θBest.m<0
            θBest.C= abs(θBest.C)
            θBest.s= abs(θBest.s)
            θBest.m= abs(θBest.m)
        end
        println("## OC fitting done.")
        return(θBest, θerr, fitfound)
    end
end

## Estimate  prior parameters
function priorGuess(rad, dens)
    m= 3
    mdisp= 2

    # King 
    # m=1.5 ; mdisp= 0.0001
    
    C= maximum(dens)
    Cdisp= C/2

    s= StatsBase.mean(rad, StatsBase.weights(convert(Array{Float64,1}, dens)))
    sdisp= s/2

    priorMean= modelCauchy(C, s, m)
    priorDisp= modelCauchy(Cdisp, sdisp, mdisp)

    return(priorMean, priorDisp)
end

####################
## radial density functions
####################

function volumeSphere(r)
    return(4 / 3 * π * r^3)
end

function areaDisk(r)
    return(π * r^2)
end

function density2D(x , y , nbin=10; norm="log")
    center= [median(x) ; median(y)]

    nxy= length(x)
    A= Array{Float64}(undef,2,nxy)
    for i in 1:nxy
        A[1,i]= x[i]
        A[2,i]= y[i]
    end

    d= Distances.Euclidean()
    r= Distances.colwise(d, A, center)

    ρ= [] ; radius= [] ; err= []
    rmax= maximum(r)
    r0=0

    if norm == "log"
        rmin= 0.5
        dlogr= (log10(rmax) - log10(rmin)) /nbin
        rlogmin= log10(rmin)
    else
        dr= rmax/nbin
    end


    for i in 1:nbin
        if norm == "log"
            r1= 10^(rlogmin+i*dlogr)
        else
            r1= i*dr
        end
        indr= (r .> r0) .& (r .<= r1)
        nstar= count(indr)
        dens= nstar/(areaDisk(r1)-areaDisk(r0))
        if nstar > 0
            errDensity= dens / sqrt(nstar)
        else
            errDensity= 1e9
        end

        push!(ρ, dens)
        push!(radius, (r0+r1)/2)
        push!(err, errDensity)
        r0= r1
    end

    return(radius, ρ , err)
end

function density3D(x , y , z, nbin=10)
    center= [median(x) ; median(y) ; median(z)]

    nxy= length(x)
    A= Array{Float64}(undef,3,nxy)
    for i in 1:nxy
        A[1,i]= x[i]
        A[2,i]= y[i]
        A[3,i]= z[i]
    end

    d= Distances.Euclidean()
    r= Distances.colwise(d, A, center)

    ρ= [] ; radius= [] ; err=[]
    rmax= maximum(r)
    dr= rmax/nbin
    r0=0
    for i in 1:nbin
        r1= i*dr
        indr= (r .> r0) .& (r .<= r1)
        nstar= count(indr)
        dens= nstar/(volumeSphere(r1)-volumeSphere(r0))
        errDensity= dens / sqrt(nstar)
        push!(ρ, dens)
        push!(radius, (r1+r0)/2)
        push!(err, errDensity)
        r0= r1
    end

    return(radius, ρ , err)
end

## Estimation of the local density around each point
## nei: number of neighbor to estimate the density
##

function locdens2d(x , y , nei=10)
    nxy= length(x)
    A= Array{Float64}(undef,2,nxy)
    center= [mean(x) ; mean(y)]

    for i in 1:nxy
        A[:,i]= [x[i] y[i]]
    end

    d= Euclidean()
    p= pairwise(d, A, A, dims=2)
    r= Distances.colwise(d, A, center)
    locdensrad= []

    for i in 1:nxy
        rsort= sort(p[i,:])
        rmax= rsort[nei+1]
        locdens= (nei+1)/areaDisk(rmax)
        push!(locdensrad, locdens)
    end

    return(r,locdensrad)
end

function locdens3d(x , y , z, nei=10)
    nxyz= length(x)
    A= Array{Float64}(undef,3,nxyz)
    center= [mean(x) ; mean(y) ; mean(z)]

    for i in 1:nxyz
        A[:,i]= [x[i] y[i] z[i]]
    end

    d= Euclidean()
    p= pairwise(d, A, A, dims=2)
    r= Distances.colwise(d, A, center)
    locdensrad= []

    for i in 1:nxyz
        rsort= sort(p[i,:])
        rmax= rsort[nei+1]
        locdens= (nei+1)/volumeSphere(rmax)
        push!(locdensrad, locdens)
    end

    return(r,locdensrad)
end


## Main MCMC to estimate parameters
## firstvalue: true(set to mean), false(random from prior)
function main_mcmc(oc::sc2dcentered, mcmc::mcmcCauchy, likelihood::Function, firstvalue=false , verbose=true)
    let
        Random.seed!()
        mci = mcCauchy(zeros(Float64,0),zeros(Float64,0),zeros(Float64,0))
        mi, probi = theta(mcmc, oc, likelihood, firstvalue)
        if verbose println("## Init done...") end

        niter = mcmc.niter
        nburn = mcmc.nburnin
        nchain = 1
        loopAgain = true
        burndone = false
        ηbest= 0 ; θbest= 0

        while loopAgain
            micurrent, probcurrent = thetaiter(mi , mcmc, oc, likelihood)

        ### Metropolis-Hasting
            α = probcurrent / probi
            if α > rand()
                mi = micurrent
                probi = probcurrent
                nchain += 1
                if (nchain%50000 == 0 && verbose ) println("### chain:",nchain) end
                if nchain > nburn && !burndone  nchain = 0 ; burndone = true end    ## burnin done
                if nchain > niter loopAgain = false end
                if burndone
                    push!(mci.C, mi.C)
                    push!(mci.s, mi.s)
                    push!(mci.m, mi.m)
                end
            else
                nchain += 1
                if (nchain%50000 == 0 && verbose) println("### chain:",nchain) end
                if nchain > nburn && !burndone   nchain = 0 ; burndone = true end    ## burnin done
                if nchain > niter loopAgain = false end
                if burndone
                    push!(mci.C, mi.C)
                    push!(mci.s, mi.s)
                    push!(mci.m, mi.m)
                end
            end

            η= likelihood(mi, oc)
            if η > ηbest
                ηbest= η
                θbest= mi
            end
        end

        if verbose println("## MCMC 2D fitting done.") end

        return(mci, θbest, ηbest)
    end
end

## density and likelihood functions
function fdens1(r, θ::modelCauchy)
    r2= r*r
    s2= θ.s*θ.s
    val = θ.C / (1+r2/s2)^θ.m
    return(val)
end

function fdens1log(r, θ::modelCauchy)
    r2= r*r
    s2= exp(θ.s)*exp(θ.s)
    val = θ.C - θ.m * log(1+r2/s2)
    return(val)
end

## using local density for all points
function likelihood2d(θ::modelCauchy, oc::sc2dcentered)
    χ2= 0
    for i in 1:oc.nxy
        desti= fdens1(oc.radius[i], θ)
        χ2+= ((desti-oc.dens[i])/  oc.densbinerr[i])^2/oc.nxy     ## σ=1
    end
    p = min(1, exp(-χ2))
    #p= min(1, 1/χ2)
    return(p)
end

### binnned 2D density
function likelihood2dbin(θ::modelCauchy, oc::sc2dcentered)
    χ2= 0
    for i in 1:oc.nbin
        desti= fdens1(oc.radbin[i], θ)
        χ2+= ((desti-oc.densbin[i]) / oc.densbinerr[i])^2
    end
    p = min(1, exp(-χ2))
    return(p)
end

### binnned 2D/3D density log
### LOG
function likelihood2dbinlog(θ::modelCauchy, oc::sc2dcentered)
    χ2= 0
    for i in 1:oc.nbin
        desti= fdens1log(oc.radbin[i], θ)
        χ2+= ((desti-oc.densbin[i]) / 1)^2
        #println(i)
        #println(θ)
        #println(oc.densbin[i])
        #println(desti)
    end
    p = min(1, exp(-χ2))

    return(p)
end


## Model iterations...

function theta(mcmc::mcmcCauchy, oc::sc2dcentered, likelihood::Function, firstvalue::Bool)
    pC  = Normal(mcmc.Cmean, mcmc.Cdisp)
    pm  = Normal(mcmc.mmean, mcmc.mdisp)
    ps  = Normal(mcmc.smean, mcmc.sdisp)

    if firstvalue
        C= mcmc.Cmean
        s= mcmc.smean
        m= mcmc.mmean
    else
        C = rand(pC)
        s = rand(ps)
        m = rand(pm)
    end

    pdfC = pdf(pC,C)
    pdfm = pdf(pm, m)
    pdfs = pdf(ps, s)

    θ= modelCauchy(C, s, m)

    pdfC = pdf(pC , C)
    pdfm = pdf(pm, m)
    pdfs = pdf(ps, s)
    ptotal = pdfC*pdfm*pdfs*likelihood(θ, oc)

    return(θ, ptotal)
end

## next iteration
function thetaiter(θi::modelCauchy, mcmc::mcmcCauchy, oc::sc2dcentered, likelihood::Function)
    pC  = Normal(mcmc.Cmean, mcmc.Cdisp)
    pm  = Normal(mcmc.mmean, mcmc.mdisp)
    ps  = Normal(mcmc.smean, mcmc.sdisp)

    C_rw= Normal(θi.C, 0.5)
    s_rw= Normal(θi.s  , 0.5)
    m_rw= Normal(θi.m , 0.5)

    new_C = 0.
    new_m = 0
    new_s = 0

    new_C  = rand(C_rw)
    new_m =  rand(m_rw)
    new_s =  rand(s_rw)

    θ= modelCauchy(new_C, new_s, new_m)

    pdfC = pdf(pC , new_C)
    pdfm = pdf(pm, new_m)
    pdfs = pdf(ps, new_s)
    ptotal = pdfC*pdfm*pdfs*likelihood(θ, oc)

    return(θ, ptotal)
end

function model_rad(radius, θ, densModel::Function)
    ρfit= []
    for i in 1:length(radius)
        ρ= densModel(radius[i], θ)
        push!(ρfit, ρ)
    end

    return(ρfit)
end
