## MCMC functions for gaia.
## 
## still in testing...
##

### generate eps, min_neightbor and min_cluster
###
### peps , pminnei , pmincl are the parameter of the prior density for the DBSCAN parameters

function theta(p::GaiaClustering.abc)
    peps         = TruncatedNormal(p.epsmean, p.epsdisp, 0.1 , 1000)
    pminnei      = TruncatedNormal(p.min_nei, p.ncoredisp, 1 , 1000.)
    pmincl      = TruncatedNormal(p.min_cl, p.ncoredisp, 1 , 1000.)
    
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
                end
            end
            if notfound 
                new_minq *= 0.9
                new_minstars = trunc(Int, 0.9 * new_minstars)
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
                    if nchain > nburn && !burndone println("### burnout done...") ; nchain = 0 ; burndone = true end
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
                    if nchain > nburn && !burndone println("### burnout done...") ; nchain = 0 ; burndone = true end
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
    println("### mini stats...")
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


function abc_mcmc_dbscan_full(dfcart::GaiaClustering.Df, params::abcfull)
    let 
        println("## ABC/MCMC for DBSCAN FULL (parameters+weighting)...")
        mci = mcfull(zeros(Float64,0),zeros(Int32,0),zeros(Int32,0) , zeros(Float64,0), zeros(Int32,0), 
        zeros(Float64,0), zeros(Float64,0), zeros(Float64,0))
        
        initial = true
        th = []
        mi = modelfull(0.0,0,0, 0., 0., 0.)
        micurrent = modelfull(0.0,0,0, 0., 0., 0.)
        probi = 0.
        probcurrent = 0.
    
        minimumQ = params.minQ
        minstars = params.minstars
        niter = params.niter
        nburn = params.nburnout
    
        println("### Minimum Q : $minimumQ")
        println("### Minimum nstars : $minstars")
        minimumQ , minstars = check_qminqstar_full(dfcart, params, minimumQ, minstars)
        println("### Minimum Q : $minimumQ")
        println("### Minimum nstars : $minstars")
        if minstars < params.forcedminstars
            minstars = params.forcedminstars
            println("### Minimum nstars forced to : $minstars") 
        end
            
        while initial
            mi, probi = theta_full(params)
            dfcartnorm = getDfcartnorm(dfcart , mi)
            qres , nstars = find_clusters(dfcartnorm, dfcart, mi)
            if qres > minimumQ && nstars >= minstars 
                println("### init done ...")
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
        end

        ministats_full(100, dfcart, mi , params)
        
        nchain = 1
        loopAgain = true
        burndone = false
                    
        while loopAgain
            micurrent, probcurrent = thetaiter_full(mi , params)
            dfcartnorm = getDfcartnorm(dfcart , micurrent)
            qres , nstars = find_clusters(dfcartnorm, dfcart, micurrent)
        
            if qres > minimumQ && nstars >= minstars 
            ### Metropolis-Hasting
                α = probcurrent / probi
                if α > rand()
                    mi = micurrent
                    probi = probcurrent
                    nchain += 1
                    if (nchain%500 == 0) println("### chain:",nchain) end
                    if nchain > nburn && !burndone println("### burnout done...") ; nchain = 0 ; burndone = true end
                    if nchain > niter loopAgain = false end
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
                    if (nchain%500 == 0) println("### chain:",nchain) end
                    if nchain > nburn && !burndone println("### burnout done...") ; nchain = 0 ; burndone = true end
                    if nchain > niter loopAgain = false end
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
        end
        println("## ABC/MCMC FULL done")
        println("##")
        return(mci)
    end
end

## mini stats over Qc and Qn
function ministats_full(niter::Int, dfcart::GaiaClustering.Df, mi::modelfull, params::abcfull)
    println("### mini stats...")
    qcmini = []
    qnmini = []
    for i in 1:niter
        mitest, probtest = thetaiter_full(mi , params)
        dfcartnorm = getDfcartnorm(dfcart , mitest)
        
        qtest , ntest = find_clusters(dfcartnorm, dfcart, mitest)
        push!(qcmini,qtest)
        push!(qnmini,ntest)
    end
    println("### Qc : ",mean(qcmini))
    println("### Qn : ",mean(qnmini))    
    
    return(0)
end


function theta_full(p::abcfull)
    peps         = TruncatedNormal(p.epsmean, p.epsdisp, 0.1 , 1000)
    pminnei      = TruncatedNormal(p.min_nei, p.ncoredisp, 1 , 1000.)
    pmincl       = TruncatedNormal(p.min_cl, p.ncoredisp, 1 , 1000.)
    pw3d         = TruncatedNormal(p.w3dmean , p.w3ddisp, 0.1, 1000)
    pwvel        = TruncatedNormal(p.wvelmean , p.wveldisp, 0.1, 1000)
    pwhrd        = TruncatedNormal(p.whrdmean , p.whrddisp, 0.1, 1000)
    
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
function thetaiter_full(θi::modelfull , p::abcfull)
    let   
        peps         = TruncatedNormal(p.epsmean, p.epsdisp, 0.1 , 1000)
        pminnei      = TruncatedNormal(p.min_nei, p.ncoredisp, 1 , 1000.)
        pmincl       = TruncatedNormal(p.min_cl, p.ncoredisp, 1 , 1000.)
        pw3d         = TruncatedNormal(p.w3dmean , p.w3ddisp, 0.1, 1000)
        pwvel        = TruncatedNormal(p.wvelmean , p.wveldisp, 0.1, 1000)
        pwhrd        = TruncatedNormal(p.whrdmean , p.whrddisp, 0.1, 1000)       
        
        eps_rw    = TruncatedNormal(θi.eps,0.4 , 0.25, 1000.)
        minnei_rw = TruncatedNormal(θi.min_nei , 4  , 1,  1000.)
        mincl_rw  = TruncatedNormal(θi.min_cl  , 4  , 1, 1000.)
        w3d_rw    = TruncatedNormal(θi.w3d,  0.4, 0.25, 1000)
        wvel_rw   = TruncatedNormal(θi.wvel, 0.4, 0.25, 1000)
        whrd_rw   = TruncatedNormal(θi.whrd, 0.4, 0.25, 1000)
        
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

#### function to adapt the minimum condition to evaluate the  solution
###
function check_qminqstar_full(dfcart::GaiaClustering.Df,
        params::abcfull, minimumQ, minstars)
    let
        new_minq = minimumQ
        new_minstars = minstars
        notfound = true
        mingoodsolution = 10   
        
        println("### Checking the minQ and minStars conditions...")
        while notfound
            goodsolutions = 0
            for i in 1:params.nburnout
                mi, probi = theta_full(params)
                dfcartnorm = getDfcartnorm(dfcart , mi)
                qres , nstars = find_clusters(dfcartnorm, dfcart, mi)
                if qres > new_minq && nstars >= new_minstars
                    goodsolutions += 1
                end
                if goodsolutions > mingoodsolution
                    notfound = false
                end
            end
            if notfound 
                new_minq *= 0.9
                new_minstars = trunc(Int, 0.9 * new_minstars)
            end
            
            if new_minstars == 0
                println("### Problem with minStars, return(0.5,5)")
                return(0.5,5)
            end      
        end
        return(new_minq, new_minstars)
    end
end