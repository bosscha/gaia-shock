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

    for cl in res.clusters
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

        if length(labels) == 0
            println("## no clusters in the DBSCAN selection...")
            return(0, 0, 0)
        end
        ### metrics of the clusters
        q2d = metric2(dfcart, labels, "spatial2d" , m.aperture2d, m.maxaperture2d, m.nboot)
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

        if m.labels == "Qc"
            println("## Selecting best cluster based on Qc")
            println("## Qc: $qc")
            println("## Qn: $nlab")
            bestlabel= findmax(qc)[2]
        elseif m.labels == "Qn"
            println("## Selecting best cluster based on Qn")
            println("## Qc: $qc")
            println("## Qn: $nlab")
            bestlabel= findmax(nlab)[2]
        elseif m.labels == "QcQn"
            println("## Selecting best cluster based on Qc and Qn")
            println("## Qc: $qc")
            println("## Qn: $nlab")
            qcqn= get_qcqn.(qc, nlab, (m , ))
            println("## QcQn: $qcqn")
            bestlabel= findmax(qcqn)[2]
        elseif m.labels == "QcQnhigh"
            println("## Selecting best cluster based on Qc and maxQn. Make sure that is what you want...")
            println("## Qc: $qc")
            println("## Qn: $nlab")
            mtemp= m
            mtemp.minQn= 30
            qcqn= get_qcqn.(qc, nlab, (mtemp , ))
            println("## QcQn: $qcqn")
            bestlabel= findmax(qcqn)[2]
        else
            println("## No rules for selecting best cluster, label 1 set by default...")
            bestlabel= 1
        end

        return(bestlabel, nlab[bestlabel], qc[bestlabel])
    end
end
## function to project (Qc,Qn) values on the best value with QcQn methd
## qn = n
## favors qc for small qn otherwise qc
function sigmoid(x, c1, c2)
    f= 1/(1+exp(-c1*(x-c2)))
    return(f)
end
function get_qcqn(qc, qn,m::GaiaClustering.meta)
    k= qc*sigmoid(qn, 0.1, m.minQn)*(1-sigmoid(qn,0.1,m.maxQn))
    # k= qc*sigmoid(qn, 0.1, m.minQn)
    return(k)
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

    debug_red(indvrad)

    if length(indvrad[indvrad]) >  1
        vrad = median(df.raw[13,indx[indvrad]])
        vraddisp = std(df.raw[13,indx[indvrad]])
    elseif length(indvrad[indvrad]) == 1
        vrad     = median(df.raw[13,indx[indvrad]])
        vraddisp = 0.
    end

    doff= sqrt(median(dfcart.data[2,indx])^2+median(dfcart.data[3,indx])^2)
    doffdeg= atand(doff/distance)

    edgeratio1, edgeratio2= edge_ratio(dfcart, indx)

    sc = SCproperties2()
    sc.nstars= nstars

    # sc.distance= 1000. / parallax

    dist, sig_dist= distance_cluster(indx, df)
    sc.distance= dist
    sc.distance_std= sig_dist
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
    sc.edgratg= edgeratio1
    sc.edgratm= edgeratio2

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

        #debugging...
        #dist_debug= median(slab[1,:])
        #debug_red("dist $dist_debug")

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
    dfsid= df.sourceid[:,setdiff(1:end,idx)]

    dfcartdata= dfcart.data[:,setdiff(1:end,idx)]
    dfcartraw= dfcart.raw[:,setdiff(1:end,idx)]
    dfcarterr= dfcart.err[:,setdiff(1:end,idx)]
    dfcartsid= dfcart.sourceid[:,setdiff(1:end,idx)]

    s=size(dfdata)

    dfnew= GaiaClustering.Df(s[2],dfdata,dfraw,dferr,dfsid)
    dfcartnew= GaiaClustering.Df(s[2],dfcartdata,dfcartraw,dfcarterr,dfcartsid)

    nrem= length(idx)
    println("### $nrem stars removed")
    return(dfnew, dfcartnew)
end

#### cycle extraction
#### namefile: prefix for plotfile
### adding boolean optim for optimization or not
function cycle_extraction_optim(df::GaiaClustering.Df, dfcart::GaiaClustering.Df, m::GaiaClustering.meta, optim=true)
    let
        println("############### cycle_extraction #########")
        println("## maximum cycle: $(m.cyclemax)")
        println("## Q metric on DBSCAN solutions: $(m.labels)")
        if optim
            println("## DBSCAN/weighting optimization. It may take time...")
        else
            println("## NO optimization, the weightings and DBSCAN parameters are fixed.")
            println("## w3d : $(m.w3d)")
            println("## wvel : $(m.wvel)")
            println("## whrd : $(m.whrd)")
            println("## ϵ : $(m.eps)")
            println("## min_cluster : $(m.mcl)")
            println("## min_neighbor : $(m.mnei)")
        end
        votname= basename(m.votname)       # get rid of the leading part only useful to read data
        cyclerun= true ; cycle= 1 ; FLAG= 0
        smallest_offdeg= 1e9 ; cycle_offdeg= -1

        sclist= [] ; mcmclist= [] ; perflist= [] ; chainlist= []

        println("##")
        while cyclerun
            FLAG= -1
            tstart= now()
            println("###############")
            print("## "); println(blue("starting cycle $cycle ..."))
            @printf("## starting time: %s \n",tstart)

            if optim
                ## extraction one cycle.. MCMC optimization
                mc , iter, FLAGmcmc= abc_mcmc_dbscan_full2(dfcart, m)
                println("## ABC/MCMC flag: $FLAGmcmc")
                nchain= length(mc.qc)
                println("## $iter iterations performed...")
                println("## $nchain chains")

                println("## optimization completed..")
                println("## analyzing solutions...")

                if length(mc.eps) > 0    ## at least 1 solution for median...
                    plot_dbscanfull_mcmc(m.plotdir, "$votname.$cycle", mc , false)

                    # __plot_mcmc_article(mc)

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
                
                    if length(labels) == 0
                        FLAGmcmc= 0   ## to force stop even w/o optimization
                        nchain= 0
                        println("### no solution found with DBSCAN...")
                    end
                else
                    FLAGmcmc= 0   ## to force stop 
                    nchain= 0
                    println("### No chain found ...")
                end
            else
                println("## setting the weightings/DBSCAN")
                eps= m.eps
                min_nei= m.mnei
                min_cl= m.mcl
                w3d= m.w3d
                wvel= m.wvel
                whrd= m.whrd

                mres = GaiaClustering.modelfull(eps,min_nei,min_cl,w3d,wvel,whrd)
                dfcartnorm = getDfcartnorm(dfcart, mres)
                labels = clusters(dfcartnorm.data ,eps  , 20, min_nei, min_cl)

                if length(labels) == 0
                    FLAGmcmc= 0   ## to force stop even w/o optimization
                    nchain= 0
                    println("### no solution found with DBSCAN...")
                else
                    FLAGmcmc= -1
                    nchain=  m.minchainreached+1
                end
            end

            if FLAGmcmc== -1 || nchain > m.minchainreached

                labelmax , nmax, qc = find_cluster_label2(labels, df, dfcart, m)

                ### tail, if step 2 extraction ("tail") is requested, here it goes...
                if m.tail == "yes" 
                    println("## Performing step 2 for extraction (tails, clumps, etc) ...")
                    println("### Cut radius: $(m.maxRadTail) (pc) -- velocity: $(m.maxVelTail) (km/s) -- CMD: $(m.maxDistCmdTail) (mag)")
                    println("### Subtracting first solution to find tails ...")

                    dfnew, dfcartnew= remove_stars(df, dfcart, labels[labelmax])
                    println("### Return the new full solution in label 1")

                    ## labelmax 1: full solution, 2: step1 solution, 3: step2 solution
                    labels, labelmax= tail_stars(df, dfcart, dfnew, dfcartnew, labels[labelmax], m, cycle=cycle)
                    nmax= size(labels[labelmax])[1]                
                end

                ## Principal components
                pc, pcres= compute_PC(df, dfcart, labels, labelmax)

                if m.pca == "no"
                    oc= export_df("$votname.$cycle", m.ocdir, df , dfcart , labels , labelmax, pc, m)
                elseif m.pca == "yes"
                    println("## PCA components added to the oc")
                    oc= export_df("$votname.$cycle", m.ocdir, df , dfcart , labels , labelmax, pc, m)
                end
                println("## label $labelmax written as an oc solution...")

                ### isochrone fitting...
                minptsiso= 10  # minimum pts to fit.
                nvalid= counting_valid(oc)
                debug_red("nvalid : $nvalid")

                if m.iso == "yes" && count(!isnan, oc.ag) > minptsiso && nvalid > minptsiso 
                    println("## Performing isochrone fitting on the solution...")
                    debug_red(count(!isnan, oc.ag))
                    oc, age, feh, feh_gaia, iso= perform_isochrone_fitting(oc, m.isomodel)
                    agemyr= 10^age / 1e6
                    total_mass= sum(oc.star_mass)

                    txt_iso= DataFrame(age=agemyr, feh= feh, feh_gaia= feh_gaia, total_mass=total_mass, cycle=cycle,  votname=votname)
                    plot_isochrone(m.plotdir, "$votname.$cycle", oc, iso, txt_iso , false)

                    name= split("$votname.$cycle",".")
                    infix= ""
                    for iname in name
                        if iname != "vot"
                            infix *= iname*"."
                        end
                    end
                    infix1 = infix * "oc.mass.csv" ; filename= @sprintf("%s/%s",m.ocdir, infix1)
                    CSV.write(filename,oc,delim=';')
                    infix2 = infix * "isochrone.csv" ; filename= @sprintf("%s/%s",m.ocdir, infix2)
                    CSV.write(filename,iso,delim=';')
                    println("## Isochrone $filename saved...")
                else
                    ### placeholder for isochrone fitting
                    agemyr= -99.9 ; feh= -99.9 ; feh_gaia= -99.9
                end
                ###

                edgeratio1, edgeratio2= edge_ratio(dfcart, labels[labelmax])
                scproperties = get_properties_SC2(labels[labelmax] , df, dfcart)
                scdf= convertStruct2Df(scproperties)
                
                ## smallest offdeg...
                if scdf.offdeg[1] < smallest_offdeg
                    smallest_offdeg= scdf.offdeg[1]
                    cycle_offdeg= cycle
                end

           
                insertcols!(scdf, 1, :votname => votname)
                insertcols!(scdf, 2, :uuid => string(m.uuid))
                insertcols!(scdf, 3, :cycle => cycle)
                insertcols!(scdf, 4, :pc3 => pcres[3])
                insertcols!(scdf, 4, :pc2 => pcres[2])
                insertcols!(scdf, 4, :pc1 => pcres[1])

                ## add solution used for DBSCAN and weighting
                insertcols!(scdf, 30, :whrd => whrd)
                insertcols!(scdf, 30, :wvel => wvel)
                insertcols!(scdf, 30, :w3d => w3d)
                insertcols!(scdf, 30, :mnei => min_nei)
                insertcols!(scdf, 30, :mcl => min_cl)
                insertcols!(scdf, 30, :eps => eps)

                ## Xg, Yg, Zg median Galactic position
                dg = df.data[3,labels[labelmax]]
                lg = df.data[1,labels[labelmax]]
                bg = df.data[2,labels[labelmax]]
                parg = df.raw[5,labels[labelmax]]
                ra = df.raw[1,labels[labelmax]]
                dec = df.raw[2,labels[labelmax]]

                xg= [] ; yg= [] ; zg= []
                for i in 1:size(dg)[1]
                    x1 = galXYZ(ra[i],dec[i],dg[i])
                    push!(xg,x1[1]) ; push!(yg,x1[2]) ; push!(zg,x1[3])
                end
                Xgm= median(xg) ; Ygm= median(yg); Zgm= median(zg)
                debug_red("Xg $Xgm Yg $Ygm Zg $Zgm")

                if m.tail == "yes"
                    ## Xg Yg Zg only for the core...
                    xgc= [] ; ygc= [] ; zgc= []
                    dgt= df.data[3,labels[2]]
                    rat = df.raw[1,labels[2]]
                    dect = df.raw[2,labels[2]]
                    for i in 1:size(dgt)[1]
                        x1 = galXYZ(rat[i],dect[i],dgt[i])
                        push!(xgc,x1[1]) ; push!(ygc,x1[2]) ; push!(zgc,x1[3])
                    end  
                    Xgm_core= median(xgc) ; Ygm_core= median(ygc); Zgm_core= median(zgc)
                    debug_red("Xg_core $Xgm_core Yg_core $Ygm_core Zg_core $Zgm_core")
                end

                ## add Galactic positions into SC
                insertcols!(scdf, 18, :Zg => Zgm)
                insertcols!(scdf, 18, :Yg => Ygm)
                insertcols!(scdf, 18, :Xg => Xgm)

                ## Galactic UVW velocities
                uvw= galUVW(scdf.ra[1], scdf.dec[1], scdf.distance[1], scdf.pmra[1], scdf.pmdec[1], scdf.vrad[1])
                insertcols!(scdf, 21, :W => uvw[3])
                insertcols!(scdf, 21, :V => uvw[2])
                insertcols!(scdf, 21, :U => uvw[1])

                ## if tail=yes add the properties of the core
                if m.tail == "yes"
                    sc_core = get_properties_SC2(labels[2] , df, dfcart)
                    ncore= size(labels[2])[1]
                    insertcols!(scdf, 33, :vraddisp_core => sc_core.vraddisp)
                    insertcols!(scdf, 33, :vbdisp_core => sc_core.vbdisp)
                    insertcols!(scdf, 33, :vldisp_core => sc_core.vldisp)                   
                    insertcols!(scdf, 33, :zdisp_core => sc_core.zdisp)
                    insertcols!(scdf, 33, :ydisp_core => sc_core.ydisp)
                    insertcols!(scdf, 33, :xdisp_core => sc_core.xdisp)
                    insertcols!(scdf, 33, :distance_core => sc_core.distance)
                    insertcols!(scdf, 33, :dec_core => sc_core.dec)
                    insertcols!(scdf, 33, :ra_core => sc_core.ra)
                    insertcols!(scdf, 33, :Zg_core => Zgm_core)
                    insertcols!(scdf, 33, :Yg_core => Ygm_core)
                    insertcols!(scdf, 33, :Xg_core => Xgm_core)
                    insertcols!(scdf, 33, :ncore => ncore)
                else
                    insertcols!(scdf, 33, :vraddisp_core => scdf.vraddisp)
                    insertcols!(scdf, 33, :vbdisp_core => scdf.vbdisp)
                    insertcols!(scdf, 33, :vldisp_core => scdf.vldisp)                   
                    insertcols!(scdf, 33, :zdisp_core => scdf.zdisp)
                    insertcols!(scdf, 33, :ydisp_core => scdf.ydisp)
                    insertcols!(scdf, 33, :xdisp_core => scdf.xdisp)
                    insertcols!(scdf, 33, :distance_core => scdf.distance)
                    insertcols!(scdf, 33, :dec_core => scdf.dec)
                    insertcols!(scdf, 33, :ra_core => scdf.ra)
                    insertcols!(scdf, 33, :ncore => scdf.nstars)
                end

                ## isochrone fitting, if not, placeholder..
                insertcols!(scdf, 7, :feh_gaia => feh_gaia)
                insertcols!(scdf, 7, :feh => feh)
                insertcols!(scdf, 7, :age => agemyr) 

                if optim
                    insertcols!(res, 2,  :cycle => cycle)
                    push!(mcmclist, res)
                    ## create DF chain
                    dfchain= create_DFchain(mc, votname, cycle)
                    push!(chainlist,dfchain)
                end

                push!(sclist, scdf)

                distance_lab= median(dg)
                distance_lab2= median(1000 ./ parg)
                debug_red("distance parallax without bstrap and weight .... $distance_lab2")

                println("###")
                println("### solution label: $labelmax")
                print("### "); println(red(@sprintf("PC1: %3.1f , PC2: %3.1f , PC3: %3.1f", pcres[1], pcres[2], pcres[3])))
                println("### Offdeg: $(scproperties.offdeg)")
                println("### Edge ratio: $(scproperties.edgratm)")
                println("### N stars: $nmax")
                println("### Distance: $distance_lab pc (no weight)")
                println("### Qc: $qc")
                println("###")

                ## score only makes sense for optimization, otherwise set to 1.
                if optim k= score_cycle(qc, nmax, nchain, iter) else k= 1. end
                @printf("## score cycle %d: %3.3f \n",cycle, k)
            
                extraplot= DataFrame(cycle=cycle, score_cycle=k, qc=qc, votname=votname, pc1=pcres[1],pc2=pcres[2], pc3=pcres[3], xg=Xgm, yg=Ygm,zg=Zgm, uuid= m.uuid)

                if m.tail == "yes"
                    labelmax= 2 #### step1 solution
                    plot_cluster2(m.plotdir, "$votname.$cycle", labels[labelmax], sc_core,
                    dfcart , false, extraplot)
                else
                    plot_cluster2(m.plotdir, "$votname.$cycle", labels[labelmax], scproperties,
                    dfcart , false, extraplot)
                end


                jump= 50  ## stars to jump in the plot
                plot_rawdata(m.plotdir, "$votname.$cycle", labels[labelmax], scproperties,
                    dfcart , pc, jump, false, extraplot)

                plot_astrom(m.plotdir, "$votname.$cycle", labels[labelmax], scproperties,
                    df , false, extraplot)

                println("###")
                println("### subtracting BEST solution from Df...")
                if m.tail == "yes"
                    println("### step2 (tail) solution is subtracted also...")
                    labelmax= 1
                end
                debug_red(size(labels[labelmax]))
                debug_red(labelmax)
                dfnew, dfcartnew= remove_stars(df, dfcart, labels[labelmax])
                
                df= dfnew
                dfcart= dfcartnew
         
                ########################### STOP conditions #########
                FLAG= 0
                if nmax < m.minstarstop
                    FLAG= FLAG | (1<<0)
                    println("### extraction stopped at cycle $cycle")
                    println("### nmax too low...")
                    cyclerun= false
                end
                if cycle == m.cyclemax
                    FLAG= FLAG | (1<<1)
                    println("### extraction stopped at cycle $cycle")
                    println("### cyclemax reached...")
                    cyclerun= false
                end
                if qc < m.qcminstop
                    FLAG= FLAG | (1<<2)
                    println("### extraction stopped at cycle $cycle")
                    println("### Qc too low...")
                    cyclerun= false
                end
                if w3d/wvel < m.wratiominstop || wvel/w3d < m.wratiominstop
                    FLAG= FLAG | (1<<3)
                    println("### extraction stopped at cycle $cycle")
                    println("### weight ratio too low...")
                    cyclerun= false
                end
                if FLAGmcmc == 3 && nchain > m.minchainreached
                    FLAG= FLAG | (1<<4)
                    println("## extraction stopped at cycle $cycle")
                    println("## chain iteration not performed completely but sufficient to keep...")
                    cyclerun= false
                end
                #####################################################
                ## Time
                ##
                tend= now()
                duration= Dates.value(tend-tstart) / (1000*1)
                nstar= size(df.raw)[2]
                if optim timeperiterstar= duration / (iter*nstar) end
                if optim timeperchainstar= duration / (nchain*nstar) end
                println("## ")
                @printf("## Time: \n")
                @printf("## duration per cycle %3.3f sec \n", duration)
                if optim @printf("## duration per iteration*star %3.3e sec \n", timeperiterstar) end
                if optim @printf("## duration per chain*star %3.3e sec \n", timeperchainstar) end
                @printf("##\n")

                if optim
                    ## log the results of performances
                    dfout= DataFrame(votname=votname, cycle=cycle, nstar=nstar, qc=qc, nmax=nmax, nchain=nchain, iter=iter,
                    scorecycle=k, duration=duration, timeperiterstar=timeperiterstar ,
                    timeperchainstar= timeperchainstar )
                    push!(perflist, dfout)
                end

                cycle += 1
            else
                println("## nothing found, stopped...")
                FLAG= 0
                cyclerun= false
            end
        end
        if cycle >= 2
            save_cycle_optim(sclist, mcmclist, perflist, chainlist, m, optim)
        end
        println("##")
        println("## Smallest offdeg: $(smallest_offdeg) for  cycle $(cycle_offdeg) ")
        println("## Copy cycle $(cycle_offdeg) to cycle 0")
        copy_cycle_0(votname,cycle_offdeg, m)

        return(cycle-1, FLAG)
    end
end
#############################
### save results cycle in csv
###
function save_cycle(sc, mcmc, perf, chain,  m::GaiaClustering.meta)
    cd(m.wdir)

    filesc= @sprintf("%s.sc.csv", m.prefile)
    filemcmc= @sprintf("%s.mcmc.csv", m.prefile)
    fileperf= @sprintf("%s.perf.csv", m.prefile)
    filechain= @sprintf("%s.chain.csv", m.prefile)
    votname= m.votname

    ncycle= length(sc)

    for i in 1:ncycle
        if !isfile(filesc)
            CSV.write(filesc,sc[i],delim=';')
            println("## $filesc created...")
            CSV.write(filemcmc,mcmc[i],delim=';')
            println("## $filemcmc created...")
            CSV.write(fileperf,perf[i],delim=';')
            println("## $fileperf created...")
            CSV.write(filechain,chain[i],delim=';')
            println("## $filechain created...")
        else
            res= CSV.File(filesc, delim=";") |> DataFrame
            append!(res,sc[i]) ; CSV.write(filesc,res,delim=';')
            res= CSV.File(filemcmc, delim=";") |> DataFrame
            append!(res,mcmc[i]) ; CSV.write(filemcmc,res,delim=';')
            res= CSV.File(fileperf, delim=";") |> DataFrame
            append!(res,perf[i]) ; CSV.write(fileperf,res,delim=';')
            # res= CSV.File(filechain, delim=";") |> DataFrame
            CSV.write(filechain,chain[i],delim=';', append=true)
        end
    end
end
##
### save results cycle in csv with the optim option
###
function save_cycle_optim(sc, mcmc, perf, chain,  m::GaiaClustering.meta,optim)
    cd(m.wdir)
    filesc= @sprintf("%s.sc.csv", m.prefile)
    if optim
        filemcmc= @sprintf("%s.mcmc.csv", m.prefile)
        fileperf= @sprintf("%s.perf.csv", m.prefile)
        filechain= @sprintf("%s.chain.csv", m.prefile)
    end
    votname= m.votname
    ncycle= length(sc)

    for i in 1:ncycle
        if !isfile(filesc)
            CSV.write(filesc,sc[i],delim=';')
            println("## $filesc created...")
            if optim
                CSV.write(filemcmc,mcmc[i],delim=';')
                println("## $filemcmc created...")
                CSV.write(fileperf,perf[i],delim=';')
                println("## $fileperf created...")
                CSV.write(filechain,chain[i],delim=';')
                println("## $filechain created...")
            end
        else
            println("## cycle $i results appended to $filesc ...")
            res= CSV.File(filesc, delim=";") |> DataFrame

            append!(res,sc[i]) ; CSV.write(filesc,res,delim=';')
            if optim
                if !isfile(filemcmc)
                    CSV.write(filemcmc,mcmc[i],delim=';')
                    println("## $filemcmc created...")
                    CSV.write(fileperf,perf[i],delim=';')
                    println("## $fileperf created...")
                    CSV.write(filechain,chain[i],delim=';')
                    println("## $filechain created...")
                else
                    res= CSV.File(filemcmc, delim=";") |> DataFrame
                    append!(res,mcmc[i]) ; CSV.write(filemcmc,res,delim=';')
                    res= CSV.File(fileperf, delim=";") |> DataFrame
                    append!(res,perf[i]) ; CSV.write(fileperf,res,delim=';')
                    # res= CSV.File(filechain, delim=";") |> DataFrame
                    CSV.write(filechain,chain[i],delim=';', append=true)
                end
            end
        end
    end
end
## compute the closeness of the cluster solution to the edge of the data.
## return a ratio in the range [0., 1]. 1 means at the edge and 0 in the center
function edge_ratio(dfcart::GaiaClustering.Df, ind)
    r2= dfcart.data[2,ind] .* dfcart.data[2,ind] .+ dfcart.data[3,ind] .* dfcart.data[3,ind]
    minX= minimum(dfcart.data[1,ind])
    maxX= maximum(dfcart.data[1,ind])
    xg= median(dfcart.data[2,ind]) ; yg= median(dfcart.data[3,ind]) ; dg= median(dfcart.data[1,ind])

    indx= (dfcart.data[1,:] .<= maxX) .& (dfcart.data[1,:] .>= minX)
    rtot2=  dfcart.data[2,indx] .* dfcart.data[2,indx]  .+ dfcart.data[3,indx] .* dfcart.data[3,indx]

    ratio= sqrt((xg^2+yg^2) / maximum(rtot2))

    alpha= maximum(sqrt.(rtot2)) / maxX
    rg= alpha*dg
    ratio_2= sqrt(maximum(r2)) / rg

    return(ratio, ratio_2)
end

## compute principal components of a salution
## compute cumulated ratio for the first 3 PCs and the 3 first PCs of the normalized data
##
function compute_PC(df::GaiaClustering.Df, dfcart::GaiaClustering.Df, labels, labelmax)
        print("### Computing principal components... \n")
        s=size(labels[labelmax])
        if s[1] == 0 
            print("### No solutions found, PCA not performed")
            return([0],[0,0,0])
        end
        data= zeros(8,s[1])

        X= dfcart.data[1, labels[labelmax]]
        Y= dfcart.data[2, labels[labelmax]]
        Z= dfcart.data[3, labels[labelmax]]
        vl= df.data[4,labels[labelmax]]
        vb= df.data[5,labels[labelmax]]
        gbar= df.data[6,labels[labelmax]]
        Crp= df.data[7,labels[labelmax]]
        Cbp= df.data[8,labels[labelmax]]

        data[1,:]= X
        data[2,:]= Y
        data[3,:]= Z
        data[4,:]= vl
        data[5,:]= vb
        data[6,:]= gbar
        data[7,:]= Crp
        data[8,:]= Cbp

        dt= StatsBase.fit(ZScoreTransform, data, dims=2)
        d2= StatsBase.transform(dt, data)
        M= fit(PCA, d2)
        p= projection(M)
        Yt = MultivariateStats.transform(M, d2)


        totvar= tvar(M)
        pvs= principalvars(M)
        ratioac= accumulate(+, pvs ./ totvar)

        if length(ratioac) >= 3
            pcres= 100 .* ratioac[1:3]
        else
            pcres= [100,100,100]
        end

        return(Yt, pcres)
end
########################################################
### copy the cycle plots with smallest offdeg to cycle 0
function copy_cycle_0(votname, cycle_offdeg, m::GaiaClustering.meta)
    debug_red("cycle 0 ...")
    debug_red(votname)
    f= glob("$(votname).$(cycle_offdeg)*png", m.plotdir)
    if size(f)[1] > 0
        for file in f
            file= splitdir(file)[2]       ### remove leading directories
            newfile= joinpath(m.plotdir,replace(file,".$(cycle_offdeg)." => ".0."))
            src= joinpath(m.plotdir, file)
            cp(src,newfile, force=true)
        end
    end
end

####################################################
### compute the distance with the bootrapped weighted parallax 
function distance_cluster(indx, df::GaiaClustering.Df)
    parallax = df.raw[5,indx]
    parallax_err = df.err[1,indx]

    dist_boot , sig_boot= bootstrapping_distance(parallax, parallax_err)

    return(dist_boot, sig_boot)
end


function bootstrapping_distance(arr, err, Nboot= 50)
    ndata= size(arr)[1]

    if Nboot > ndata
        Nboot= trunc(Int, ndata / 2)
        println("### Warning Nbootstrapping down to $Nboot ..")
    end

    nsample= 100000
    idx= 1:ndata
    res= zeros(nsample)

    println("## Bootstrapping the star parallax with $nsample samples of $Nboot...")

    for i in 1:nsample
        indsh= shuffle(idx)
        para_s= arr[indsh[1:Nboot]]
        para_err_s= err[indsh[1:Nboot]]

        weight= 1 ./ (para_err_s .* para_err_s)
        norm= sum(weight)
        weight = weight ./ norm
        res[i]= sum(weight .* para_s)
    end

    dist= median(1000 ./ res)
    sig= std(1000 ./ res)
    println("### Distance: $dist (std: $sig) pc")


    return(dist, sig)
end