## MCMC functions for gaia.
##
##

### generate eps, min_neightbor and min_cluster
###
### peps , pminnei , pmincl are the parameter of the prior density for the DBSCAN parameters

function theta(p::GaiaClustering.abc)
    peps         = TruncatedNormal(p.epsmean, p.epsdisp, 0.1 , 1000)
    pminnei      = TruncatedNormal(p.min_nei, p.ncoredisp, 1 , 1000.)
    pmincl       = TruncatedNormal(p.min_cl, p.ncoredisp, 1 , 1000.)

    e = rand(peps)          ; pe = pdf(peps,e)
    n = trunc(Int,rand(pminnei))  ; pn = pdf(pminnei, n)
    c = trunc(Int,rand(pmincl))   ; pc = pdf(pmincl, c)

    ptotal = pe*pn*pc
    params = GaiaClustering.model(e, n, c)
    return(params, ptotal)
end

## iterate with random walk and yield the probability
##
function thetaiter(θi::GaiaClustering.model , p::GaiaClustering.abc)
    let
        peps         = TruncatedNormal(p.epsmean, p.epsdisp, 0.1 , 1000)
        pminnei      = TruncatedNormal(p.min_nei, p.ncoredisp, 1 , 1000.)
        pmincl       = TruncatedNormal(p.min_cl, p.ncoredisp, 1 , 1000.)

        ## Random Walk settings..
        eps_rw    = TruncatedNormal(θi.eps,0.25 , 0.1, 1000.)
        minnei_rw = TruncatedNormal(θi.min_nei , 4  , 1,  1000.)
        mincl_rw  = TruncatedNormal(θi.min_cl  , 4  , 1, 1000.)

        new_e = 0.
        new_mn = 0
        new_mcl = 0

        iternotfound =true
        while iternotfound
            new_e = rand(eps_rw)
            new_mn =  trunc(Int,rand(minnei_rw))
            new_mcl = trunc(Int,rand(mincl_rw))

            if new_e > 0 && new_mn > 0 && new_mcl > 0
                iternotfound = false
            end
        end

        pe = pdf(peps , new_e)
        pn = pdf(pminnei, new_mn)
        pc = pdf(pmincl, new_mcl)
        ptotal = pe * pn * pc

        params = GaiaClustering.model(new_e, new_mn, new_mcl)
        return(params, ptotal)
    end
end

#### function to adapt the minimum condition to evaluate the  solution
###
function check_qminqstar(df::GaiaClustering.Df, dfcart::GaiaClustering.Df,
        params::GaiaClustering.abc, minimumQ, minstars)
    let
        new_minq = minimumQ
        new_minstars = minstars
        notfound = true
        mingoodsolution = 50

        println("### Checking the minQ and minStars conditions...")
        println("### Minimum good solutions $mingoodsolution")
        while notfound
            goodsolutions = 0
            for i in 1:params.nburnout
                mi, probi = theta(params)
                qres , nstars = find_clusters(df, dfcart, mi)
                if qres > new_minq && nstars >= new_minstars
                    goodsolutions += 1
                end
                if goodsolutions > mingoodsolution
                    notfound = false
                    return(new_minq, new_minstars)
                end
            end
            if notfound
                new_minq *= 0.9
                new_minstars = trunc(Int, 0.9 * new_minstars)
                println("### MinQ not reached yet... testing with $new_minq")
            end

            if new_minstars == 0
                println("### Problem with minStars, return(0.5,5)")
                return(0.5,5)
            end
        end
        return(new_minq, new_minstars)
    end
end

## ABC MCMC (following Weyan et al. 2013)
## first simple scheme for testing
##
function abc_mcmc_dbscan(df::GaiaClustering.Df, dfcart::GaiaClustering.Df, params::GaiaClustering.abc)
    let
        println("## ABC/MCMC for DBSCAN...")
        mci = GaiaClustering.mc(zeros(Float64,0),zeros(Int32,0),zeros(Int32,0) , zeros(Float64,0), zeros(Int32,0))

        initial = true
        th = []
        mi = GaiaClustering.model(0.0,0,0)
        micurrent = GaiaClustering.model(0.0,0,0)
        probi = 0.
        probcurrent = 0.

        minimumQ = params.minQ
        minstars = params.minstars
        niter = params.niter
        nburn = params.nburnout

        println("### Minimum Q : $minimumQ")
        println("### Minimum nstars : $minstars")
        minimumQ , minstars = check_qminqstar(df, dfcart, params, minimumQ, minstars)
        #params.minQ = minimuQ
        #params.minstars = minstars
        println("### Minimum Q : $minimumQ")
        println("### Minimum nstars : $minstars")
        println("### Lower conditions lead to lower quality solutions...")
        println("###")

        while initial
            mi, probi = theta(params)
            qres , nstars = find_clusters(df, dfcart, mi)
            if qres > minimumQ && nstars >= minstars
                println("### init done ...")
                initial = false
                push!(mci.eps, mi.eps)
                push!(mci.mne, mi.min_nei)
                push!(mci.mcl, mi.min_cl)
                push!(mci.qc,  qres)
                push!(mci.qn,  nstars)
            end
        end

        ministats(100, df, dfcart, mi , params)

        nchain = 1
        loopAgain = true
        burndone = false

        while loopAgain
            micurrent, probcurrent = thetaiter(mi , params)
            qres , nstars = find_clusters(df, dfcart, micurrent)

            if qres > minimumQ && nstars >= minstars
            ### Metropolis-Hasting
                α = probcurrent / probi
                if α > rand()
                    mi = micurrent
                    probi = probcurrent
                    nchain += 1
                    if (nchain%500 == 0) println("### chain:",nchain) end
                    if nchain > nburn && !burndone println("### burn-in done...") ; nchain = 0 ; burndone = true end
                    if nchain > niter loopAgain = false end
                    if burndone
                        push!(mci.eps, mi.eps)
                        push!(mci.mne, mi.min_nei)
                        push!(mci.mcl, mi.min_cl)
                        push!(mci.qc , qres)
                        push!(mci.qn,  nstars)
                    end
                else
                    nchain += 1
                    if (nchain%500 == 0) println("### chain:",nchain) end
                    if nchain > nburn && !burndone println("### burn-in done...") ; nchain = 0 ; burndone = true end
                    if nchain > niter loopAgain = false end
                    if burndone
                        push!(mci.eps, mi.eps)
                        push!(mci.mne, mi.min_nei)
                        push!(mci.mcl, mi.min_cl)
                        push!(mci.qc , qres)
                        push!(mci.qn,  nstars)
                    end
                end
            end
        end
        println("## ABC/MCMC done")
        println("##")
        return(mci)
    end
end

## mini stats over Qc and Qn
function ministats(niter::Int, df::GaiaClustering.Df, dfcart::GaiaClustering.Df,mi::GaiaClustering.model, params::GaiaClustering.abc)
    println("\n### mini stats...")
    qcmini = []
    qnmini = []
    for i in 1:niter
        mitest, probtest = thetaiter(mi , params)
        qtest , ntest = find_clusters(df, dfcart, mitest)
        push!(qcmini,qtest)
        push!(qnmini,ntest)
    end
    println("### Qc : ",mean(qcmini))
    println("### Qn : ",mean(qnmini))

    return(0)
end

##############################################
#### MCMC FULL for dbscan
######
function getDfcartnorm(dfcart::GaiaClustering.Df, mc::modelfull)
    blck       = [[1,2,3],[4,5], [6,7,8]]
    wghtblck   = [mc.w3d, mc.wvel, mc.whrd]
    norm       = "identity"

    dfcartnorm , scale8 = normalization_PerBlock(dfcart, blck, wghtblck , norm, false, false)
    return(dfcartnorm)
end

## optimize dbscan parameters with ABC/MCMC method
##
function abc_mcmc_dbscan_full2(dfcart::GaiaClustering.Df, params::GaiaClustering.meta)
    let
        Random.seed!()
        println("## ABC/MCMC for DBSCAN FULL (parameters+weighting)...")
        println("## ABC/MCMC v2")

        mci = mcfull(zeros(Float64,0),zeros(Int32,0),zeros(Int32,0) , zeros(Float64,0), zeros(Int32,0),
        zeros(Float64,0), zeros(Float64,0), zeros(Float64,0))

        initial = true
        FLAG= -1            ## everything normal
        th = []
        mi = modelfull(0.0,0,0, 0., 0., 0.)
        micurrent = modelfull(0.0,0,0, 0., 0., 0.)
        probi = 0.
        probcurrent = 0.

        minimumQ= params.minQc
        minstars= params.minQn
        maxstars= params.maxQn
        niter= params.nchain             ## number of chains
        nburn= params.nburnout

        maxiter= params.maxiter          ## twice, for init and normal iteration

        println("### Burn-in : $nburn")
        println("### Chains  : $niter")
        println("### Minimum Qc : $minimumQ")
        println("### Minimum Qn : $minstars")
        println("### Maximum Qn : $maxstars")
        println("### Maximum iterations: $maxiter")
        println("###")

        minimumQ , minstars = check_qminqstar_full2(dfcart, params)
        @printf("### Minimum Qc: %3.3f \n",minimumQ)
        println("### Minimum Qn: $minstars")
        if minstars < params.forcedminstars
            minstars = params.forcedminstars
            println("### Minimum Qn forced to $minstars")
        end

        iter= 0
        p = Progress(maxiter,"Initialization...")

        while initial
            mi, probi = theta_full(params)
            dfcartnorm = getDfcartnorm(dfcart , mi)
            qres , nstars = find_clusters2(dfcartnorm, dfcart, mi, params)
            if qres > minimumQ && nstars >= minstars && nstars <= maxstars
                println("\n### init done ...")
                initial = false
                push!(mci.eps, mi.eps)
                push!(mci.mne, mi.min_nei)
                push!(mci.mcl, mi.min_cl)
                push!(mci.w3d, mi.w3d)
                push!(mci.wvel, mi.wvel)
                push!(mci.whrd, mi.whrd)
                push!(mci.qc,  qres)
                push!(mci.qn,  nstars)
            end
            iter += 1
            update!(p, iter)
            if iter > maxiter
                println("\n### Maximum iteration reached, no solution for init...")
                FLAG= 2      ## init not performed completely
                return(mci, iter,FLAG)
            end
        end

        ministats_full(100, dfcart, mi , params)

        nchain= 1
        loopAgain= true
        burndone= false
        tstart= now()

        iter= 0
        p = Progress(niter, "Optimizing DBSCAN for Extraction...")

        while loopAgain
            micurrent, probcurrent = thetaiter_full(mi , params)
            dfcartnorm = getDfcartnorm(dfcart , micurrent)
            qres , nstars = find_clusters2(dfcartnorm, dfcart, micurrent, params)

            if qres > minimumQ && nstars >= minstars && nstars <= maxstars
            ### Metropolis-Hasting
                α = probcurrent / probi
                if α > rand()
                    mi = micurrent
                    probi = probcurrent
                    nchain += 1
                    if (nchain%1000 == 0) println("\n### chain:",nchain) end
                    if nchain > nburn && !burndone println("\n### burn-in done...") ; nchain = 0 ; burndone = true end
                    if nchain > niter && burndone loopAgain = false end
                    if burndone
                        push!(mci.eps, mi.eps)
                        push!(mci.mne, mi.min_nei)
                        push!(mci.mcl, mi.min_cl)
                        push!(mci.w3d, mi.w3d)
                        push!(mci.wvel, mi.wvel)
                        push!(mci.whrd, mi.whrd)
                        push!(mci.qc , qres)
                        push!(mci.qn,  nstars)
                    end
                else
                    nchain += 1
                    if (nchain%1000 == 0) println("\n### chain:",nchain) end
                    if nchain > nburn && !burndone println("\n### burn-in done...") ; nchain = 0 ; burndone = true end
                    if nchain > niter && burndone loopAgain = false end
                    if burndone
                        push!(mci.eps, mi.eps)
                        push!(mci.mne, mi.min_nei)
                        push!(mci.mcl, mi.min_cl)
                        push!(mci.w3d, mi.w3d)
                        push!(mci.wvel, mi.wvel)
                        push!(mci.whrd, mi.whrd)
                        push!(mci.qc , qres)
                        push!(mci.qn,  nstars)
                    end
                end
            end

            iter += 1
            update!(p, nchain)
            if (iter%10000 == 0)
                titer= now()
                duration= Dates.value(titer-tstart) / (1000*3600)
                meanTime= duration / nchain
                eta= @sprintf("%3.3f",meanTime * (niter-nchain))
                println("\n### iteration: $iter (ETA:$eta h)")
            end
            if iter > maxiter
                println("\n### Maximum iteration reached, current solution returned...")
                if burndone   FLAG= 3 end    ## nchain not reached before maxiter but burnin done...
                if !burndone  FLAG= 4 end    ## burn-in not performed
                return(mci, iter, FLAG)
            end
        end
        println("## ABC/MCMC (full) done")
        println("##")
        return(mci,  iter, FLAG)
    end
end

## mini stats over Qc and Qn
function ministats_full(niter::Int, dfcart::GaiaClustering.Df, mi::modelfull, params::meta)
    println("\n### mini stats...")
    qcmini = []
    qnmini = []
    for i in 1:niter
        mitest, probtest = thetaiter_full(mi , params)
        dfcartnorm = getDfcartnorm(dfcart , mitest)

        qtest , ntest = find_clusters2(dfcartnorm, dfcart, mitest, params)
        push!(qcmini,qtest)
        push!(qnmini,ntest)
    end
    @printf("### Qc : %3.3f \n",mean(qcmini))
    @printf("### Qn : %3.3f \n",mean(qnmini))

    debug_red("ministats_full")
    return(0)
end


function theta_full(p::meta)
    peps         = truncated(Normal(p.epsmean, p.epsdisp), 0.1 , 1000)
    pminnei      = truncated(Normal(p.min_nei, p.ncoredisp), 1 , 1000.)
    pmincl       = truncated(Normal(p.min_cl, p.ncoredisp), 1 , 1000.)
    pw3d         = truncated(Normal(p.w3dmean , p.w3ddisp), 0.1, 1000)
    pwvel        = truncated(Normal(p.wvelmean , p.wveldisp), 0.1, 1000)
    pwhrd        = truncated(Normal(p.whrdmean , p.whrddisp), 0.1, 1000)

    e = rand(peps)          ; pe = pdf(peps,e)
    n = trunc(Int,rand(pminnei))  ; pn = pdf(pminnei, n)
    c = trunc(Int,rand(pmincl))   ; pc = pdf(pmincl, c)
    w3 = rand(pw3d)          ; pw3 = pdf(pw3d,w3)
    wv = rand(pwvel)         ; pwv = pdf(pwvel,wv)
    wh = rand(pwhrd)         ; pwh = pdf(pwhrd,wh)

    ptotal = pe*pn*pc*pw3*pwv*pwh
    params = modelfull(e, n, c , w3, wv, wh)
    return(params, ptotal)
end


## iterate with random walk and yield the probability
##
function thetaiter_full(θi::modelfull , p::meta)
    let
        peps         = truncated(Normal(p.epsmean, p.epsdisp), 0.1, 1000)
        pminnei      = truncated(Normal(p.min_nei, p.ncoredisp), 1 , 1000.)
        pmincl       = truncated(Normal(p.min_cl, p.ncoredisp), 1 , 1000.)
        pw3d         = truncated(Normal(p.w3dmean , p.w3ddisp), 0.1, 1000)
        pwvel        = truncated(Normal(p.wvelmean , p.wveldisp), 0.1, 1000)
        pwhrd        = truncated(Normal(p.whrdmean , p.whrddisp), 0.1, 1000)

        eps_rw    = truncated(Normal(θi.eps,p.eps_rw) , 0.1, 1000.)
        minnei_rw = truncated(Normal(θi.min_nei , p.minnei_rw) , 1,  1000.)
        mincl_rw  = truncated(Normal(θi.min_cl  , p.mincl_rw)  , 1, 1000.)
        w3d_rw    = truncated(Normal(θi.w3d,  p.w3d_rw), 0.1, 1000)
        wvel_rw   = truncated(Normal(θi.wvel, p.wvel_rw), 0.1, 1000)
        whrd_rw   = truncated(Normal(θi.whrd, p.whrd_rw), 0.1, 1000)

        new_e = 0.
        new_mn = 0
        new_mcl = 0
        new_w3d = 0.
        new_wvel = 0.
        new_whrd = 0.

        iternotfound = true
        while iternotfound
            new_e    = rand(eps_rw)
            new_mn   = trunc(Int,rand(minnei_rw))
            new_mcl  = trunc(Int,rand(mincl_rw))
            new_w3d  = rand(w3d_rw)
            new_wvel = rand(wvel_rw)
            new_whrd = rand(whrd_rw)

            ### adding the condition on minclusterand > minnei..
            if new_e > 0 && new_mn > 0 && new_mcl > 0 && new_mn <= new_mcl
                iternotfound = false
            end
        end

        pe = pdf(peps , new_e)
        pn = pdf(pminnei, new_mn)
        pc = pdf(pmincl, new_mcl)
        pw3 = pdf(pw3d, new_w3d)
        pwv = pdf(pwvel, new_wvel)
        pwh = pdf(pwhrd, new_whrd)

        ptotal = pe * pn * pc * pw3 * pwv * pwh
        params = modelfull(new_e, new_mn, new_mcl , new_w3d, new_wvel, new_whrd)
        return(params, ptotal)
    end
end

## check maximum iterations
function check_qminqstar_full2(dfcart::GaiaClustering.Df, params::GaiaClustering.meta)
    let
        new_minq = params.minQc
        new_minstars = params.minQn
        notfound = true

        ## standard sofar 50
        mingoodsolution = params.mingoodsolution
        niter = params.niterqminq
        maxiter= niter*30           ## cycle numbers * niter

        println("#### Checking the minQc and minQn conditions...")
        println("#### Minimum good solutions: $mingoodsolution")
        println("#### Number of iterations: $niter, maxiter: $maxiter")

        totaliter= 0
        while notfound
            goodsolutions = 0
            for i in 1:niter
                # println("-- $i")
                mi, probi = theta_full(params)
                # println(mi)
                # println(probi)
                dfcartnorm = getDfcartnorm(dfcart , mi)
                qres , nstars = find_clusters2(dfcartnorm, dfcart, mi, params)
                if qres > new_minq && nstars >= new_minstars
                    goodsolutions += 1
                end
                if goodsolutions > mingoodsolution
                    notfound = false
                    return(new_minq, new_minstars)
                end
            end
            totaliter += niter
            if notfound
                new_minq *= 0.95
                new_minstars = trunc(Int, 0.95 * new_minstars)
                new_minqstr= @sprintf("%3.3f", new_minq)
                println("#### MinQ not reached yet... testing with $new_minqstr")
            end

            if new_minstars == 0 || totaliter > maxiter
                println("#### Maximum iteration reached..,")
                if new_minstars < 5
                    new_minstars= 5
                end
                if new_minq < 0.5
                    new_minq= 0.5
                end
                println("#### Maximum iterations reached.. minimums set to $new_minq, $new_minstars")
                return(new_minq, new_minstars)
            end
        end
        return(new_minq, new_minstars)
    end
end

function extraction_mcmc(mc, votname)
    epsm = median(mc.eps)
    epsd = std(mc.eps)
    mneim = median(mc.mne)
    mneid = std(mc.mne)
    mclm = median(mc.mcl)
    mcld = std(mc.mcl)
    qcm = median(mc.qc)
    qnm = median(mc.qn)
    qcd = std(mc.qc)
    qnd = std(mc.qn)
    w3dm = median(mc.w3d)
    w3dd = std(mc.w3d)
    wvelm = median(mc.wvel)
    wveld = std(mc.wvel)
    whrdm = median(mc.whrd)
    whrdd = std(mc.whrd)

    @printf("## DBSCAN/MCMC stats: \n")
    @printf("### ϵ : %3.3f +/- %3.3f \n", epsm, epsd)
    @printf("### min_nei  : %3.1f +/- %3.3f \n", mneim, mneid)
    @printf("### min_clus  : %3.1f +/- %3.3f \n", mclm, mcld)
    @printf("### W3d  : %3.3f +/- %3.3f \n", w3dm, w3dd)
    @printf("### Wvel  : %3.3f +/- %3.3f \n" , wvelm, wveld)
    @printf("### Whrd  : %3.3f +/- %3.3f \n", whrdm, whrdd)
    @printf("### Qn  : %3.3f +/- %3.3f \n",qnm, qnd)
    @printf("### Qc  : %3.3f +/- %3.3f \n",qcm, qcd)
    @printf("##")

    res = DataFrame(votname=votname, epsm = epsm, epsd=epsd, mneim=mneim,mneid=mneid,mclm=mclm,mcld=mcld,
            qcm=qcm,qcd=qcd, qnm=qnm,qnd=qnd,
            w3dm=w3dm,w3dd=w3dd,wvelm=wvelm,wveld=wveld,whrdm=whrdm,whrdd=whrdd)
    return(res)
end

### create DataFrame for the chaine
###
function create_DFchain(mc, votname, cycle)
    df= DataFrame(votname=votname, cycle=cycle, eps=mc.eps , mne=mc.mne, mcl=mc.mcl, w3d=mc.w3d, wvel=mc.wvel, whrd=mc.whrd, qc=mc.qc, qn=mc.qn)
    return(df)
end
