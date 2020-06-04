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
function find_clusters(df::GaiaClustering.Df, dfcart::GaiaClustering.Df , m::GaiaClustering.model,
    aperture2d = 1.5, maxaperture2d = 15, aperturev = 3.0, maxaperturev = 20, nboot = 50 ,
    aperture3d = 3., maxaperture3d = 20)
    let
        println("Warning...find_clusters deprecated, find_clusters2 should be used instead")
        labels = clusters(df.data , m.eps , 20, m.min_nei, m.min_cl)
        nsol= length(labels)
        if nsol == 0
            return(0, 0)
        end
        if nsol > 1000
            println("### Warning... $nsol clusters found in find_clusters ")
        end

    ### metrics of the clusters
        q2d = metric2(dfcart, labels, "spatial2d" , aperture2d, maxaperture2d, nboot)
        q3d = metric2(dfcart, labels, "spatial3d" , aperture3d, maxaperture3d, nboot)     #### Added
        qv  = metric2(dfcart, labels, "velocity" , aperturev, maxaperturev, nboot)
        qp, qa = metric2(dfcart, labels, "HRD" )

        nlab = []
        for ilab in labels
            push!(nlab,length(ilab))
        end


    #### metric for the number of stars in the cluster
        qn = []
        for nl in nlab
            push!(qn,log10(nl))
        end

        qc = 0.
        qstar = 0
        for i in 1:length(nlab)
            k1 = q2d[i][1]
            k1bis = q3d[i][1]
            k2 = qv[i][1]
            k3 = qa[i][1]
            k4 = qn[i]
            ############### Composite metric ###
            qq = (2k1 + k1bis + 3k2 + k3 + k4) / 8.0
            # qq = (3k1 + k1bis + 3k2 + k4) / 8.0
            ###############
            if qq > qc
                qc = qq
                qstar = nlab[i]
            end
        end
        return(qc, qstar)
    end
end

## broadcasting modelfull ..
##
function find_clusters(df::GaiaClustering.Df, dfcart::GaiaClustering.Df , m::GaiaClustering.modelfull,
    aperture2d = 1.5, maxaperture2d = 15, aperturev = 3.0, maxaperturev = 20, nboot = 50,
    aperture3d = 3., maxaperture3d = 20)
    let
        println("Warning...find_clusters deprecated, find_clusters2 should be used instead")
        labels = clusters(df.data , m.eps , 20, m.min_nei, m.min_cl)
        nsol= length(labels)
        if nsol == 0 || nsol > 1000
            println("### Warning... $nsol clusters found in find_clusters2 ")
            return(0, 0)
        end

    ### metrics of the clusters
        q2d = metric2(dfcart, labels, "spatial2d" , aperture2d, maxaperture2d, nboot)
        q3d = metric2(dfcart, labels, "spatial3d" , aperture3d, maxaperture3d, nboot)
        qv  = metric2(dfcart, labels, "velocity" , aperturev, maxaperturev, nboot)
        qp, qa = metric2(dfcart, labels, "HRD" )

        nlab = []
        for ilab in labels
            push!(nlab,length(ilab))
        end


    #### metric for the number of stars in the cluster
        qn = []
        for nl in nlab
            push!(qn,log10(nl))
        end

        qc = 0.
        qstar = 0
        for i in 1:length(nlab)
            k1 = q2d[i][1]
            k1bis = q3d[i][1]
            k2 = qv[i][1]
            k3 = qa[i][1]
            k4 = qn[i]
            ############### Composite metric ###
            qq = (2k1 + k1bis + 3k2 + k3 + k4) / 8.0
            # qq = (3k1 + k1bis + 3k2 + k4) / 8.0
            ###############
            if qq > qc
                qc = qq
                qstar = nlab[i]
            end
        end
        return(qc, qstar)
    end
end
