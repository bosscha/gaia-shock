## function for the 2nd step of extraction.
## Apply a cut in radius, velocity and a match with the first CMD 
## dfnew: df removed from the step1 solution
## testing stars with very similar cmd
function tail_stars(df::GaiaClustering.Df, dfcart::GaiaClustering.Df, dfnew::GaiaClustering.Df, dfcartnew::GaiaClustering.Df, idx,  
    m::GaiaClustering.meta ; cycle=1, plot=true)
    debug_red("entering  tail...")
    doc= transform_df(df, dfcart, idx)
    dfstep1= doc   ## step1 before filtering
    
    ## source id for the full solution, starting with the old one...

    s=dfnew.ndata
    s1= length(dfnew.data[4,:])
    dnew= transform_df(dfnew, dfcartnew, 1:s)

    debug_red(median(dnew.Y))

    xcenter= median(doc.X) ; ycenter= median(doc.Y) ;  zcenter= median(doc.Z)
    vlc= median(doc.vl) ; vbc= median(doc.vb) ; vrc= median(doc.vrad)

    radiusMax= m.maxRadTail  ## distance max to oc
    r2= radiusMax*radiusMax

    velocityMax= m.maxVelTail ## velocity difference max
    v2= velocityMax*velocityMax

    dnewrad=filter(row -> ((row.X-xcenter)*(row.X-xcenter)+(row.Y-ycenter)*(row.Y-ycenter)+(row.Z-zcenter)*(row.Z-zcenter)) < r2, dnew)
    dnewvel=filter(row -> ((row.vl .- vlc)*(row.vl .- vlc) + (row.vb .- vbc)*(row.vb .- vbc)) < v2, dnewrad)
    
    s= length(dfnew.data[1, :])
    srad= length(dnewrad.X)
    svel= length(dnewvel.X)

    debug_red("Filtering..")
    debug_red("$s $srad $svel")

    doc= filter(:BmR0 => x -> !(ismissing(x) || isnothing(x) || isnan(x)), doc)
    debug_red("doc filtered $(size(doc))")
    dnewvel= filter(:BmR0 => x -> !(ismissing(x) || isnothing(x) || isnan(x)), dnewvel)
    
    idx_tail , dist= distance_cmd_tail(doc,dnewvel)

    # __plot_dist_cmd(dist)
    
    ## cut with cmd
    cmdDistMax= m.maxDistCmdTail
    idc= findall(x->(x< cmdDistMax),dist)

    dnewcmd= dnewvel[idc,:]

    sid_final= vcat(df.sourceid[1,idx],dnewcmd.sourceid) 
    
    label_newsolution= []

    sold= size(idx)[1] ; snew= size(idc)[1]
    println("### Tail solutions step 1: $sold")
    println("### Tail solutions step 2: $snew")

    for sid in sid_final
        inew= findall(x-> x == sid, df.sourceid[1,:])
        if length(inew) > 1 || length(inew) == 0
            println("### Warning,there is duplicated id  in the solution ... ")
        else
            push!(label_newsolution, inew[1])
        end
    end
    
    ## labels for step 2 solution
    label_step2= []
    for sid in dnewcmd.sourceid
        inew= findall(x-> x == sid, df.sourceid[1,:])
        if length(inew) > 1 || length(inew) == 0
            println("### Warning,there is duplicated id  in the solution ... ")
        else
            push!(label_step2, inew[1])
        end
    end

    ## final results
    dfres= transform_df(df, dfcart,label_newsolution)
   
    dfstep2= dnewcmd
    dfstep2= transform_df(df, dfcart,label_step2)

    labels= [label_newsolution, idx, label_step2]     ## index for solution
    labelmax= 1                                       ## solution is label 1, step 1 -> label 2, step 2 -> label 3

    if plot 
        println("### Plot the full results...")
        votname= basename(m.votname)
        nstep1= sold[1] ; nstep2= snew[1] ; ntotal= nstep1+nstep2
        
        ## fit density2D to Cauchy
        nbin= 25
        fit, err, found=  spatialParameter("", nbin=nbin, verbose=false, niter=15000, dfoc=dfres)    
        fit1, err1, found1=  spatialParameter("", nbin=nbin, verbose=false, niter=15000, dfoc=dfstep1)

        dfinfo= DataFrame(cycle=cycle, nstep1=nstep1, nstep2=nstep2, ntotal=ntotal)
        pc=[]
        oc= export_df("$votname.$cycle", m.ocdir, df , dfcart , labels , labelmax, pc, m, save=false)

        ## plot_tail should be cleaned!! 
        plot_tail(m.plotdir, votname, dfres, dfstep1, dfstep2 ,  dist, fit, err, found, fit1, err1, found1, dfinfo, oc)
    end
   
    return(labels,  labelmax)
end

function transform_df(df,dfcart, idx)
    X= dfcart.data[1, idx]
    Y= dfcart.data[2, idx]
    Z= dfcart.data[3, idx]
    dist= df.data[3,idx]
    vl= df.data[4,idx]
    vb= df.data[5,idx]
    vrad= df.raw[13,idx]
    gbar= df.raw[10,idx]
    rp= df.raw[11,idx]
    bp= df.raw[12,idx]
    ag= df.raw[14,idx]
    a0= df.raw[15,idx]
    ebmr= df.raw[16,idx]
    source_id= df.sourceid[1,idx]

    gbar0 = gbar .- ag
    BmR0 = bp .- rp .- ebmr
    ## absolute magnitude
    G = gbar0 .+ 5 .- 5log10.(dist)

    dft= DataFrame(X=X,Y=Y,Z=Z,dist=dist,vl=vl,vb=vb,vrad=vrad,gbar=gbar,rp=rp,bp=bp,ag=ag,a0=a0,ebmr=ebmr,gbar0=gbar0,BmR0=BmR0, G=G, sourceid=source_id)

    return(dft)
end

##### distance of two CMDs
#### df1, df2: two solutions
#### cut= maximum distance to cut
function distance_cmd_tail(df1, df2)
    n= length(df1.X)
    dfref= zeros(2,n)
    dfref[1,:] = df1.BmR0
    dfref[2,:] = df1.G

    kdtree = KDTree(dfref; leafsize = 20)

    debug_red("knn dref: $n")

    npts= length(df2.G)
    pts= zeros(2,npts)
    pts[1,:] = df2.BmR0
    pts[2,:] = df2.G

    if n>1
        idxs, dists = nn(kdtree, pts)
        idx= collect(Iterators.flatten(idxs))
        dist= collect(Iterators.flatten(dists))
        idx=collect(1:npts)
    else
        println("## Warning: Reference tree is too small for KNN analysis of the CMD distance, return a fake one...")
        idx= [1,2,3,4]
        dist=[1e6,1e6,1e6,1e6]
    end

    return(idx , dist)
end

function _density_count(xx, yy, nbin=128, xrange=[-100,100],yrange=[-100,100])
    data = (xx,yy) 
    stepx= (xrange[2]-xrange[1])/nbin
    stepy= (yrange[2]-yrange[1])/nbin 
    debug_red("density..")

    h = FHist.Hist2D((xx,yy), (xrange[1]:stepx:xrange[2], yrange[1]:stepy:yrange[2]))
    dens= h.hist.weights ./ (stepx*stepy)

    wav= atrous(dens, 7)
    # wavFilt= thresholdingWav(wav,Normal())
    rec= addWav(wav,4,8)
    rec= permutedims(rec, [2, 1])
    nrec= size(rec)

    nticks= 4
    xti= [] ; yti= [] ; xv=[] ; yv=[]
    dx= (xrange[2]-xrange[1])/nticks  ; dy= (yrange[2]-yrange[1])/nticks
    dxi= (nrec[1]-1)/nticks ; dyi= (nrec[2]-1)/nticks

    for i in 1:(nticks+1)
        xx= xrange[1] + (i-1)*dx  ; sx= @sprintf("%3.0f",xx) ; xi= (i-1)*dxi
        yy= yrange[1] + (i-1)*dy  ; sy= @sprintf("%3.0f",yy) ; yi= (i-1)*dyi
        push!(xti, sx) ; push!(xv, xi)
        push!(yti, sy) ;  push!(yv, yi) 
    end

    PyPlot.plt.figure(figsize=(9.0,8.0))
    ax= PyPlot.plt.subplot(1, 1, 1 )

    ax.set_xticks(xv) ;  ax.set_xticklabels(xti)
    ax.set_yticks(yv) ;  ax.set_yticklabels(yti)

    vmin, vmax=  __level_dens(rec , 2, 5)
    nlev= 20
    PyPlot.plt.contour(rec, nlev, vmin=vmin, vmax=vmax, linewidths= 0.2, colors= "black") 

    PyPlot.plt.savefig("test_density.png")
    return(dens)
end
