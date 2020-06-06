### functions to analyze stellar clusters
###
### HR diagram metric using Delaunay tessalation to estimate density in the HRD
###
## The tesselation algorithm works with coordinates in range [1,2]
#
function metricHRD(hrdi::Df , label)

    QA = []
    QP = []

    for ilab in label
        nhr = size(ilab)
        points = []

        for n in 1:nhr[1]
            ## Fixed 27.02.2020!!
            push!(points,[hrdi.data[7,ilab[n]],-hrdi.data[6,ilab[n]]])
        end

        ## println("Voronoi tesselation ...")
        ## min_pts is the minimum of pts
        if length(points) > 4
            per, area = voronoi(points)
        else
            per  = [1e10,1e10,1e10,1e10]
            area = [1e10,1e10,1e10,1e10]
        end

        p = sort(per ./ nhr[1])
        a = sort(area ./ nhr[1])

        n10 = max(1, trunc(Int, nhr[1] / 10))
        start10 = max(1 , 1*n10)
        end10   = min(nhr[1] , 9*n10)

        q1a  = mean(a[start10:end10])
        q1p  = mean(p[start10:end10])
        dq1a = std(a[start10:end10])
        dq1p = std(p[start10:end10])

        qa  = -log(max(1e-15, q1a))
        dqa = abs(dq1a/qa)
        qp  = -log(max(1e-15, q1p))
        dqp = abs(dq1p / qp)

        push!(QA, (qa,dqa))
        push!(QP, (qp,dqp))

    end
    return(QA,QP)
end


##  proj :
##  spatial2d : y,z (id 2,3)
##  spatial3d : x,y,z (id 1,2,3)
##  velocity  : vrad,vdec (id 4,5)
##  HRD       : Diagram HR (id 6,7)
##
## labels : labels of the clusters##
##
## !!! works with dfcart !!! (cartesian and normalized)

function metric(s::Df, labels ,proj = "spatial2d", APERTURE = 1.0 , MAXAPERTURE = 15.0, NBOOTSTRAP = 50)

    if proj == "spatial2d"
        ycenter = mean(s.data[2,:])
        zcenter = mean(s.data[3,:])
    elseif proj == "spatial3d"
        xcenter = mean(s.data[1,:])
        ycenter = mean(s.data[2,:])
        zcenter = mean(s.data[3,:])
    elseif proj == "velocity"
        vxcenter = mean(s.data[4,:])
        vycenter = mean(s.data[5,:])
    elseif proj == "HRD"
        res = metricHRD(s, labels)
        return(res)
    end

    aper2 = APERTURE * APERTURE
    Q = []
    for ilab in labels
        slab = s.data[:,ilab]

        if proj == "spatial2d"
            yy = slab[2,:]
            zz = slab[3,:]
            ycenter = mean(yy)
            zcenter = mean(zz)
            angle_out = 2π * rand(Float64, NBOOTSTRAP)
            rad_out   = (MAXAPERTURE - 2APERTURE) * rand(Float64, NBOOTSTRAP)
            angle_in  = 2π * rand(Float64, NBOOTSTRAP)
            rad_in    = APERTURE * randexp(Float64, NBOOTSTRAP)
        elseif proj == "spatial3d"
            xx = slab[1,:]
            yy = slab[2,:]
            zz = slab[3,:]
            xcenter = mean(xx)
            ycenter = mean(yy)
            zcenter = mean(zz)
            angle_out = 2π * rand(Float64, NBOOTSTRAP)
            phi_out   = π * rand(Float64, NBOOTSTRAP)
            rad_out   = (MAXAPERTURE - 2APERTURE) * rand(Float64, NBOOTSTRAP)
            angle_in  = 2π * rand(Float64, NBOOTSTRAP)
            phi_in    = π * rand(Float64, NBOOTSTRAP)
            rad_in    = APERTURE * randexp(Float64, NBOOTSTRAP)
        elseif proj == "velocity"
            vx = slab[4,:]
            vy = slab[5,:]
            vxcenter = mean(vx)
            vycenter = mean(vy)
            angle_out = 2π * rand(Float64, NBOOTSTRAP)
            rad_out   = (MAXAPERTURE - 2APERTURE) * rand(Float64, NBOOTSTRAP)
            angle_in  = 2π * rand(Float64, NBOOTSTRAP)
            rad_in    = APERTURE * rand(Float64, NBOOTSTRAP)
        end

        qc = zeros(NBOOTSTRAP)

        for k in 1:NBOOTSTRAP
            if proj == "spatial2d"
                yout = ycenter + rad_out[k] * cos(angle_out[k])
                zout = zcenter + rad_out[k] * sin(angle_out[k])
                radii_out = (slab[2,:] .- yout) .* (slab[2,:] .- yout) .+ (slab[3,:] .- zout) .* (slab[3,:] .- zout)
                yin = ycenter + rad_in[k] * cos(angle_in[k])
                zin = zcenter + rad_in[k] * sin(angle_in[k])
                radii_in = (slab[2,:] .- yin) .* (slab[2,:] .- yin) .+ (slab[3,:] .- zin) .* (slab[3,:] .- zin)
            elseif proj == "spatial3d"
                xout = xcenter + rad_out[k] * cos(angle_out[k]) *  cos(phi_out[k])
                yout = ycenter + rad_out[k] * sin(angle_out[k]) *  cos(phi_out[k])
                zout = zcenter + rad_out[k] * sin(phi_out[k])
                radii_out = (slab[1,:] .- xout) .* (slab[1,:] .- xout) .+ (slab[2,:] .- yout) .* (slab[2,:] .- yout).+ (slab[3,:] .- zout) .* (slab[3,:] .- zout)
                xin = xcenter + rad_in[k] * cos(angle_in[k]) *  cos(phi_in[k])
                yin = ycenter + rad_in[k] * sin(angle_in[k]) *  cos(phi_in[k])
                zin = zcenter + rad_in[k] * sin(phi_in[k])
                radii_in = (slab[1,:] .- xin) .* (slab[1,:] .- xin) .+ (slab[2,:] .- yin) .* (slab[2,:] .- yin).+ (slab[3,:] .- zin) .* (slab[3,:] .- zin)
            elseif proj == "velocity"
                vxout = vxcenter + rad_out[k] * cos(angle_out[k])
                vyout = vycenter + rad_out[k] * sin(angle_out[k])
                radii_out = (slab[4,:] .- vxout) .* (slab[4,:] .- vxout) .+ (slab[5,:] .- vyout) .* (slab[5,:] .- vyout)
                vxin = vxcenter + rad_in[k] * cos(angle_in[k])
                vyin = vycenter + rad_in[k] * sin(angle_in[k])
                radii_in = (slab[4,:] .- vxin) .* (slab[4,:] .- vxin) .+ (slab[5,:] .- vyin) .* (slab[5,:] .- vyin)
            end

            isdout = radii_out .< aper2
            nout   = length(radii_out[isdout])
            isin   = radii_in .< aper2
            nin    = length(radii_in[isin])

            qc[k] = log(max(1,nin)) - log(max(1,nout))

        end

        Qc    = mean(qc)
        Q_std = std(qc)
        #println("Q: $Qc -- Q_std : $Q_std")
        push!(Q,(Qc,Q_std))
    end
    return(Q)
end

### DBSCAN
function clusters(data , epsilon, leaf , minneigh, mincluster)
    eps = epsilon
    leafsize = leaf
    min_neighbors = minneigh
    min_cluster_size = mincluster

    res = dbscan(data , eps , leafsize = leaf, min_neighbors = minneigh, min_cluster_size=mincluster)

    label = Vector{Vector{Int}}()

    for cl in res
        indx = cl.core_indices
        append!(indx, cl.boundary_indices)
        push!(label,indx)
    end
    return(label)
end

## find clusters from dbscan with a metric
##
## compute the dbscan clusters, the metric and the Metropolis acceptance
##
## df: Gaia struct (data.jl), normally a weighted cartersian df (dfcartnorm)
## dfcart: Gaia struct in cartesian w/o any transformation (dfcart)
##
## !! *qq* : composite metric. The weight could be adjusted.
## !! *qstar* : number of stars in the highest qq
## m: dbscan parameters.
##
## WARNING: the default values for APERTURE
## Check the metric(2) function used...
## Note 2.3.2020: the HRD metric is discarded and the weight scaling is down 6


### find_clusters with meta parameters.
###
function find_clusters2(df::GaiaClustering.Df, dfcart::GaiaClustering.Df , m::GaiaClustering.modelfull,
    param::GaiaClustering.meta, verbose= false)
    let
        labels = clusters(df.data , m.eps , 20, m.min_nei, m.min_cl)
        nsol= length(labels)
        if nsol == 0 || nsol > param.clustermax
            if verbose println("### Warning... $nsol clusters found with find_clusters2, returns 0 ") end
            return(0, 0)
        end


    ### metrics of the clusters
        q2d = metric2(dfcart, labels, "spatial2d" , param.aperture2d, param.maxaperture2d, param.nboot)
        q3d = metric2(dfcart, labels, "spatial3d" , param.aperture3d, param.maxaperture3d, param.nboot)
        qv  = metric2(dfcart, labels, "velocity" , param.aperturev, param.maxaperturev, param.nboot)
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
## label from the dbscan labels with maximum stars
function find_cluster_label(labels)
    let
    println("## Selecting best cluster based on Nmax...")
    i = 1 ; nmax = 0 ; ilabel = 0
        for ilab in labels
            nlab = length(ilab)
            if nlab > nmax
                ilabel = i
                nmax= nlab
            end
            i += 1
        end
    return(ilabel, nmax)
    end
end

## label with higher Qc is selected
function find_cluster_label2(labels, df::GaiaClustering.Df, dfcart::GaiaClustering.Df ,
    m::GaiaClustering.meta)
    let
        println("## Selecting best cluster based on Qc..")
        ### metrics of the clusters
        q2d = metric2(dfcart, labels, "spatial2d" , m.aperture2d, m.maxaperture2d, m,nboot)
        q3d = metric2(dfcart, labels, "spatial3d" , m.aperture3d, m.maxaperture3d, m.nboot)     #### Added
        qv  = metric2(dfcart, labels, "velocity" , m.aperturev, m.maxaperturev, m.nboot)
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

        qc= []
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
            push!(qc,qq)
        end

        println("## Qc: $qc")
        bestlabel= findmax(qc)[2]

        return(bestlabel, nlab[bestlabel], qc[bestlabel])
    end
end

## compute the properties of the cluster with indices indx
function get_properties_SC(indx, df::GaiaClustering.Df, dfcart::GaiaClustering.Df)::SCproperties
    nstars   = length(indx)
    distance = mean(df.data[3,indx])
    l        = mean(df.data[1,indx])
    b        = mean(df.data[2,indx])
    ra       = mean(df.raw[1,indx])
    dec      = mean(df.raw[2,indx])
    vl       = mean(df.data[4,indx])
    vb       = mean(df.data[5,indx])
    xdisp    = std(dfcart.data[1,indx])
    ydisp    = std(dfcart.data[2,indx])
    zdisp    = std(dfcart.data[3,indx])
    vldisp   = std(df.data[4,indx])
    vbdisp   = std(df.data[5,indx])
    parallax = mean(df.raw[5,indx])
    pmra     = mean(df.raw[6,indx])
    pmdec    = mean(df.raw[7,indx])
    pml      = mean(df.raw[8,indx])
    pmb      = mean(df.raw[9,indx])

    vrad     = 0.
    vraddisp = 0.
    indvrad = isnotnan(df.raw[13,indx])

    if length(indvrad) >  1
        vrad = mean(df.raw[13,indx[indvrad]])
        vraddisp = std(df.raw[13,indx[indvrad]])
    elseif length(indvrad) == 1
        vrad     = mean(df.raw[13,indx[indvrad]])
        vraddisp = 0.
    end

    sc = SCproperties(nstars , distance, ra , dec , l , b , parallax, pmra , pmdec , pml, pmb, vl , vb, vrad, xdisp ,
        ydisp , zdisp , vldisp, vbdisp, vraddisp)
    return(sc)

end

## compute the properties of the cluster with indices indx
## new extended version of SCproperties.
## mean -> median
function get_properties_SC2(indx, df::GaiaClustering.Df, dfcart::GaiaClustering.Df)::SCproperties2
    nstars   = length(indx)
    distance = median(df.data[3,indx])
    l        = median(df.data[1,indx])
    b        = median(df.data[2,indx])
    ra       = median(df.raw[1,indx])
    dec      = median(df.raw[2,indx])
    vl       = median(df.data[4,indx])
    vb       = median(df.data[5,indx])
    xdisp    = std(dfcart.data[1,indx])
    ydisp    = std(dfcart.data[2,indx])
    zdisp    = std(dfcart.data[3,indx])
    vldisp   = std(df.data[4,indx])
    vbdisp   = std(df.data[5,indx])
    parallax = median(df.raw[5,indx])
    pmra     = median(df.raw[6,indx])
    pmdec    = median(df.raw[7,indx])
    pml      = median(df.raw[8,indx])
    pmb      = median(df.raw[9,indx])

    vrad     = 0.
    vraddisp = 0.
    indvrad = isnotnan(df.raw[13,indx])

    if length(indvrad) >  1
        vrad = mean(df.raw[13,indx[indvrad]])
        vraddisp = std(df.raw[13,indx[indvrad]])
    elseif length(indvrad) == 1
        vrad     = mean(df.raw[13,indx[indvrad]])
        vraddisp = 0.
    end

    doff= sqrt(median(dfcart.data[2,indx])^2+median(dfcart.data[3,indx])^2)
    doffdeg= atand(doff/distance)

    sc = SCproperties2()

    sc.nstars= nstars
    sc.distance= distance
    sc.ra= ra
    sc.dec= dec
    sc.l= l
    sc.b= b
    sc.parallax= parallax
    sc.pmra= pmra
    sc.pmdec= pmdec
    sc.pml= pml
    sc.pmb= pmb
    sc.vl= vl
    sc.vb= vb
    sc.vrad= vrad
    sc.xdisp= xdisp
    sc.ydisp= ydisp
    sc.zdisp= zdisp
    sc.vldisp= vldisp
    sc.vbdisp= vbdisp
    sc.vraddisp= vraddisp
    sc.offdeg= doffdeg

    return(sc)
end
## some of the parameters are changed in find_clusters..
function metric2(s::GaiaClustering.Df, labels ,proj = "spatial2d", APERTURE = 1.0 ,
        MAXAPERTURE = 15.0, NBOOTSTRAP = 50 , MAXDISP2D = 10.,  MAXDISP3D = 30. , MAXDISPVEL = 4.)

    ### probabilities of the dispersion..
    p2d  = Normal(0., MAXDISP2D)
    p3d  = Normal(0., MAXDISP3D)
    pvel = Normal(0., MAXDISPVEL)

    if proj == "spatial2d"
        ycenter = mean(s.data[2,:])
        zcenter = mean(s.data[3,:])
    elseif proj == "spatial3d"
        xcenter = mean(s.data[1,:])
        ycenter = mean(s.data[2,:])
        zcenter = mean(s.data[3,:])
    elseif proj == "velocity"
        vxcenter = mean(s.data[4,:])
        vycenter = mean(s.data[5,:])
    elseif proj == "HRD"
        res = metricHRD(s, labels)
        return(res)
    end

    aper2 = APERTURE * APERTURE
    Q = []
    for ilab in labels
        slab = s.data[:,ilab]

        if proj == "spatial2d"
            yy = slab[2,:]
            zz = slab[3,:]
            ycenter = mean(yy)
            zcenter = mean(zz)
            angle_out = 2π * rand(Float64, NBOOTSTRAP)
            rad_out   = (MAXAPERTURE - 2APERTURE) * rand(Float64, NBOOTSTRAP)
            angle_in  = 2π * rand(Float64, NBOOTSTRAP)
            rad_in    = APERTURE * randexp(Float64, NBOOTSTRAP)
            ydisp = std(slab[2,:])
            zdisp = std(slab[3,:])
            prob2d = pdf(p2d,sqrt(ydisp^2+zdisp^2)) / pdf(p2d,0.)
        elseif proj == "spatial3d"
            xx = slab[1,:]
            yy = slab[2,:]
            zz = slab[3,:]
            xcenter = mean(xx)
            ycenter = mean(yy)
            zcenter = mean(zz)
            angle_out = 2π * rand(Float64, NBOOTSTRAP)
            phi_out   = π * rand(Float64, NBOOTSTRAP)
            rad_out   = (MAXAPERTURE - 2APERTURE) * rand(Float64, NBOOTSTRAP)
            angle_in  = 2π * rand(Float64, NBOOTSTRAP)
            phi_in    = π * rand(Float64, NBOOTSTRAP)
            rad_in    = APERTURE * randexp(Float64, NBOOTSTRAP)
            xdisp = std(slab[1,:])
            ydisp = std(slab[2,:])
            zdisp = std(slab[3,:])
            prob3d = pdf(p3d,sqrt(xdisp^2+ydisp^2+zdisp^2)) / pdf(p3d,0.)
        elseif proj == "velocity"
            vx = slab[4,:]
            vy = slab[5,:]
            vxcenter = mean(vx)
            vycenter = mean(vy)
            angle_out = 2π * rand(Float64, NBOOTSTRAP)
            rad_out   = (MAXAPERTURE - 2APERTURE) * rand(Float64, NBOOTSTRAP)
            angle_in  = 2π * rand(Float64, NBOOTSTRAP)
            rad_in    = APERTURE * rand(Float64, NBOOTSTRAP)
            vxdisp = std(slab[4,:])
            vydisp = std(slab[5,:])
            probvel = pdf(pvel,sqrt(vxdisp^2+vydisp^2)) / pdf(pvel,0.)
        end

        qc = zeros(NBOOTSTRAP)

        for k in 1:NBOOTSTRAP
            if proj == "spatial2d"
                prob =  prob2d
                yout = ycenter + rad_out[k] * cos(angle_out[k])
                zout = zcenter + rad_out[k] * sin(angle_out[k])
                radii_out = (slab[2,:] .- yout) .* (slab[2,:] .- yout) .+ (slab[3,:] .- zout) .* (slab[3,:] .- zout)
                yin = ycenter + rad_in[k] * cos(angle_in[k])
                zin = zcenter + rad_in[k] * sin(angle_in[k])
                radii_in = (slab[2,:] .- yin) .* (slab[2,:] .- yin) .+ (slab[3,:] .- zin) .* (slab[3,:] .- zin)
            elseif proj == "spatial3d"
                prob = prob3d
                xout = xcenter + rad_out[k] * cos(angle_out[k]) *  cos(phi_out[k])
                yout = ycenter + rad_out[k] * sin(angle_out[k]) *  cos(phi_out[k])
                zout = zcenter + rad_out[k] * sin(phi_out[k])
                radii_out = (slab[1,:] .- xout) .* (slab[1,:] .- xout) .+ (slab[2,:] .- yout) .* (slab[2,:] .- yout).+ (slab[3,:] .- zout) .* (slab[3,:] .- zout)
                xin = xcenter + rad_in[k] * cos(angle_in[k]) *  cos(phi_in[k])
                yin = ycenter + rad_in[k] * sin(angle_in[k]) *  cos(phi_in[k])
                zin = zcenter + rad_in[k] * sin(phi_in[k])
                radii_in = (slab[1,:] .- xin) .* (slab[1,:] .- xin) .+ (slab[2,:] .- yin) .* (slab[2,:] .- yin).+ (slab[3,:] .- zin) .* (slab[3,:] .- zin)
            elseif proj == "velocity"
                prob = probvel
                vxout = vxcenter + rad_out[k] * cos(angle_out[k])
                vyout = vycenter + rad_out[k] * sin(angle_out[k])
                radii_out = (slab[4,:] .- vxout) .* (slab[4,:] .- vxout) .+ (slab[5,:] .- vyout) .* (slab[5,:] .- vyout)
                vxin = vxcenter + rad_in[k] * cos(angle_in[k])
                vyin = vycenter + rad_in[k] * sin(angle_in[k])
                radii_in = (slab[4,:] .- vxin) .* (slab[4,:] .- vxin) .+ (slab[5,:] .- vyin) .* (slab[5,:] .- vyin)
            end

            isdout = radii_out .< aper2
            nout   = length(radii_out[isdout])
            isin   = radii_in .< aper2
            nin    = length(radii_in[isin])

            qc[k] = prob * (log(max(1,nin)) - log(max(1,nout)))

        end

        Qc    = mean(qc)
        Q_std = std(qc)
        #println("Q: $Qc -- Q_std : $Q_std")
        push!(Q,(Qc,Q_std))
    end
    return(Q)
end


### get the metrics
### still a list of label lists...
##
function get_metrics(labels, dfcart::GaiaClustering.Df , m::meta)
    let
        #### metrics of the clusters
        q2d = metric2(dfcart, labels, "spatial2d" , m.aperture2d, m.maxaperture2d, m.nboot)
        q3d = metric2(dfcart, labels, "spatial3d" , m.aperture3d, m.maxaperture3d, m.nboot)
        qv  = metric2(dfcart, labels, "velocity" ,  m.aperturev, m.maxaperturev, m.nboot)
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

#### cycle extraction
#### namefile: prefix for plotfile
function cycle_extraction(df::GaiaClustering.Df, dfcart::GaiaClustering.Df, m::GaiaClustering.meta)
    let
        println("############### cycle_extraction #########")

        votname= m.votname
        cyclerun= true ; cycle= 1 ; FLAG= 0

        # cyclemax= m.cyclemax
        # minstarselection=   m.minstarselection    # minimum of stars to select solution in a cycle...????
        # minstarstop=   m.minstarstop         # condition to stop cycling
        # minchainreached=  m.minchainreached      # minimum chain to analyze solution
        # qcmin=  m.qcmin                # more condition on Qc to stop cycling after the first
        # wratiomin=  m.wratiomin          # minimum ratio btwn w3d and wvel (otherwise not an OC)

        println("##")
        while cyclerun
            FLAG= -1
            tstart= now()
            println("####################")
            println("## starting cycle $cycle ...")
            @printf("## starting time: %s \n",tstart)
            ## extraction one cycle.. MCMC optimization
            mc , iter, FLAGmcmc= abc_mcmc_dbscan_full2(dfcart, m)
            println("## ABC/MCMC flag: $flag")
            nchain= length(mc.qc)
            println("## $iter iterations performed...")
            println("## $nchain chains")

            if FLAGmcmc== -1 || nchain > m.minchainreached
                println("## optimization completed..")
                println("## analyzing solutions...")
                plot_dbscanfull_mcmc(m.plotdir, "$votname.$cycle", mc , false)

                ## get the cluster and plot it
                println("## extracting the cluster using DBSCAN/WEIGHTING with:")
                res= extraction_mcmc(mc, m.votname)
                eps= res.epsm[1]
                min_nei= trunc(Int,res.mneim[1] + 0.5)
                min_cl= trunc(Int,res.mclm[1] + 0.5)
                w3d= res.w3dm[1]
                wvel= res.wvelm[1]
                whrd= res.whrdm[1]

                mres = GaiaClustering.modelfull(eps,min_nei,min_cl,w3d,wvel,whrd)
                dfcartnorm = getDfcartnorm(dfcart, mres)
                labels = clusters(dfcartnorm.data ,eps  , 20, min_nei, min_cl)
                labelmax , nmax, qc = find_cluster_label2(labels, df, dfcart, m)
                println("## label $labelmax written to oc...")
                export_df("$votname.$cycle", m.ocdir, df , dfcart, labels , labelmax)

                scproperties = get_properties_SC2(labels[labelmax] , df, dfcart)
                plot_cluster2(plotdir, "$votname.$cycle", labels[labelmax], scproperties,  dfcart , false)

                println("###")
                println("### label solution: $labelmax")
                println("### N stars: $nmax")
                println("### Qc: $qc")
                println("###")

                k= score_cycle(qc, nmax, nchain, iter)
                @printf("## score cycle %d: %3.3f \n",cycle, k)

                println("###")
                println("### subtracting BEST solution from Df...")
                dfnew, dfcartnew= remove_stars(df, dfcart, labels[labelmax])
                df= dfnew
                dfcart= dfcartnew

                ########################### STOP conditions #########
                FLAG= 0
                if nmax < m.minstarstop
                    FLAG += 1<<1
                    println("### extraction stopped at cycle $cycle")
                    println("### nmax too low...")
                    cyclerun= false
                end
                if cycle == m.cyclemax
                    FLAG += 1<<2
                    println("### extraction stopped at cycle $cycle")
                    println("### cyclemax reached...")
                    cyclerun= false
                end
                if qc < m.qcmin
                    FLAG += 1<<3
                    println("### extraction stopped at cycle $cycle")
                    println("### Qc too low...")
                    cyclerun= false
                end
                if w3d/wvel < m.wratiomin || wvel/w3d < m.wratiomin
                    FLAG += 1<<4
                    println("### extraction stopped at cycle $cycle")
                    println("### weight ratio too low...")
                    cyclerun= false
                end
                if FLAGmcmc == 3 && nchain > m.minchainreached
                    FLAG += 1<<5
                    println("## extraction stopped at cycle $cycle")
                    println("## chain iteration not performed completely but sufficient to keep...")
                    cyclerun= false
                end
                #####################################################
                ## Time
                tend= now()
                duration= Dates.value(tend-tstart) / (1000*1)
                nstar= size(df.raw)[2]
                timeperiterstar= duration / (iter*nstar)
                timeperchainstar= duration / (nchain*nstar)
                @printf("## Time: \n")
                @printf("## duration per cycle %3.3f sec \n", duration)
                @printf("## duration per iteration*star %3.3e sec \n", timeperiterstar)
                @printf("## duration per chain*star %3.3e sec \n", timeperchainstar)
                @printf("##\n")

                ## log the results of performances
                dfout= DataFrame(votname=votname, cycle=cycle, nstar=nstar, qc=qc, nmax=nmax, nchain=nchain, iter=iter,
                scorecycle=k, duration=duration, timeperiterstar=timeperiterstar ,
                timeperchainstar= timeperchainstar )
                # _updt!(filedebug, dfout)
                cycle += 1
            else
                println("## nothing found, stopped...")
                FLAG= 0
                cyclerun= false
            end
        end
        return(cycle-1, FLAG)
    end
end

function score_cycle(qc, qn, nchain, iter)
    k= log10(qc*qn*nchain /iter)
    return(k)
end

function remove_stars(df::GaiaClustering.Df, dfcart::GaiaClustering.Df, idx)
    s0=size(df.data)

    diff= setdiff(1:s0[2],idx)
    dfdata= df.data[:,setdiff(1:end,idx)]
    dfraw= df.raw[:,setdiff(1:end,idx)]
    dferr= df.err[:,setdiff(1:end,idx)]

    dfcartdata= dfcart.data[:,setdiff(1:end,idx)]
    dfcartraw= dfcart.raw[:,setdiff(1:end,idx)]
    dfcarterr= dfcart.err[:,setdiff(1:end,idx)]

    s=size(dfdata)

    dfnew= GaiaClustering.Df(s[2],dfdata,dfraw,dferr)
    dfcartnew= GaiaClustering.Df(s[2],dfcartdata,dfcartraw,dfcarterr)

    nrem= length(idx)
    println("### $nrem stars removed")
    return(dfnew, dfcartnew)
end
