## deprecated functions..

function abc_mcmc_dbscan_full(dfcart::GaiaClustering.Df, params::abcfull)
    let
        println("# abc_mcmc_dbscan_full is deprecated ...")

        Random.seed!()
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

        maxstars= 5000

        println("### Minimum Q : $minimumQ")
        println("### Minimum nstars : $minstars")
        println("### Maximum nstars : $maxstars")
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
            if qres > minimumQ && nstars >= minstars && nstars <= maxstars
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

            if qres > minimumQ && nstars >= minstars && nstars <= maxstars
            ### Metropolis-Hasting
                α = probcurrent / probi
                if α > rand()
                    mi = micurrent
                    probi = probcurrent
                    nchain += 1
                    if (nchain%1000 == 0) println("### chain:",nchain) end
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
                    if (nchain%1000 == 0) println("### chain:",nchain) end
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

### function to adapt the minimum condition to evaluate the  solution
###
function check_qminqstar_full(dfcart::GaiaClustering.Df,
        params::abcfull, minimumQ, minstars)
    let
        new_minq = minimumQ
        new_minstars = minstars
        notfound = true

        ## standard sofar 50
        mingoodsolution = 10
        niter = 500

        println("#### Checking the minQ and minStars conditions...")
        println("#### Minimum good solutions $mingoodsolution")

        while notfound
            goodsolutions = 0
            for i in 1:niter
                mi, probi = theta_full(params)
                dfcartnorm = getDfcartnorm(dfcart , mi)
                qres , nstars = find_clusters(dfcartnorm, dfcart, mi)
                if qres > new_minq && nstars >= new_minstars
                    goodsolutions += 1
                end
                if goodsolutions > mingoodsolution
                    notfound = false
                    return(new_minq, new_minstars)
                end
            end
            if notfound
                new_minq *= 0.95
                new_minstars = trunc(Int, 0.95 * new_minstars)
                println("#### MinQ not reached yet... testing with $new_minq")
            end

            if new_minstars == 0
                println("#### Problem with minStars, return(0.5,5)")
                return(0.5,5)
            end
        end
        return(new_minq, new_minstars)
    end
end
