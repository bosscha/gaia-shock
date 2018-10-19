## MCMC functions for gaia.
## 
## still in testing...
##

### generate eps, min_neightbor and min_cluster
###
### peps , pminnei , pmincl are the parameter of the prior density for the DBSCAN parameters

function _theta(p::GaiaClustering.abc)
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
function _thetaiter(θi::GaiaClustering.model , p::GaiaClustering.abc)
    let   
        peps         = TruncatedNormal(p.epsmean, p.epsdisp, 0.1 , 1000)
        pminnei      = TruncatedNormal(p.min_nei, p.ncoredisp, 1 , 1000.)
        pmincl       = TruncatedNormal(p.min_cl, p.ncoredisp, 1 , 1000.)
        
        eps_rw    = TruncatedNormal(θi.eps,0.5 , 0.1, 1000.)
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
        mingoodsolution = 10   ## !!!! Check out!!!!
        
        println("### Checking the minQ and minStars conditions...")
        while notfound
            goodsolutions = 0
            for i in 1:params.nburnout
                mi, probi = _theta(params)
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
            mi, probi = _theta(params)
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
            micurrent, probcurrent = _thetaiter(mi , params)
            qres , nstars = find_clusters(df, dfcart, micurrent)
        
            if qres > minimumQ && nstars >= minstars 
            ### Metropolis-Hasting
                α = probcurrent / probi
                if α > rand() 
                    nchain += 1
                    if (nchain%500 == 0) println("### chain:",nchain) end
                    if nchain > nburn && !burndone println("### burnout done...") ; nchain = 0 ; burndone = true end
                    if nchain > niter loopAgain = false end
                    if burndone
                        mi = micurrent
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
        mitest, probtest = _thetaiter(mi , params)
        qtest , ntest = find_clusters(df, dfcart, mitest)
        push!(qcmini,qtest)
        push!(qnmini,ntest)
    end
    println("### Qc : ",mean(qcmini))
    println("### Qn : ",mean(qnmini))    
    
    return(0)
end