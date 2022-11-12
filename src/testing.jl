## testing function (for notebooks e.g.)

function __plot_check(dfcart,plotdir,plotfile, showplot=true)
    cart= DataFrame(X=dfcart.data[1,:], Y=dfcart.data[2,:], Z=dfcart.data[3,:])

    println("## check plot subtraction ...")

    PyPlot.plt.figure(figsize=(9.0,8.0))
    PyPlot.plt.subplot(1, 1, 1 , xlim=[-100,100])
    PyPlot.plt.scatter(cart.Y, cart.X, s = 0.1 )
    PyPlot.plt.xlabel("Y (pc)")
    PyPlot.plt.ylabel("X (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.savefig(plotdir*"/"*plotfile)
    if showplot PyPlot.plt.show() end
end

function __plot_nstars(nstarh,plotfile="test-stats-votable.png", plotdir= ".")
    println("## plotting distribution...")
    PyPlot.plt.figure(figsize=(9.0,8.0))
    PyPlot.plt.subplot(1, 1, 1 )
    nbins = 50
    PyPlot.plt.hist(nstarh,nbins, range = [0,3e5],  color = "g", alpha=0.8 , label = "Votable stars",density=false)
    PyPlot.plt.xlabel("Stars")
    PyPlot.plt.ylabel("N")
    PyPlot.plt.grid(true)
    PyPlot.plt.savefig(plotdir*"/"*plotfile)
    PyPlot.plt.show()
end

function __plot_tail(df, doc, filename)
    PyPlot.plt.figure(figsize=(12.0,10.0))

    PyPlot.plt.subplot(2, 2, 1 , xlim=[-100,100])
    PyPlot.plt.scatter(df.Y, df.X, s = 0.1 )
    PyPlot.plt.xlabel("Y (pc)")
    PyPlot.plt.ylabel("X (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(2, 2, 2 , xlim=[-100,100])
    PyPlot.plt.scatter(df.Z, df.X, s = 0.1 )
    PyPlot.plt.xlabel("Z (pc)")
    PyPlot.plt.ylabel("X (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(2, 2, 3 , xlim=[-100,100] , ylim=[-100,100])
    PyPlot.plt.scatter(df.Y, df.Z, s = 0.1 )
    PyPlot.plt.xlabel("Y (pc)")
    PyPlot.plt.ylabel("Z (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(2, 2, 4 , ylim= [12, -2] )
    PyPlot.plt.scatter(df.BmR0, df.G, s = 0.1 )
    PyPlot.plt.scatter(doc.BmR0, doc.G, s = 0.1 , color="r")
    PyPlot.plt.xlabel("B-R")
    PyPlot.plt.ylabel("G")
    PyPlot.plt.grid(true)

    PyPlot.plt.savefig(filename)
end


function __plot_dist_cmd(dist, plotfile="test-dist_cmd.png", plotdir= ".")
    PyPlot.plt.figure(figsize=(9.0,8.0))
    PyPlot.plt.subplot(1, 1, 1 )
    nbins = 50
    PyPlot.plt.hist(dist,nbins, range = [0,0.50],  color = "g", alpha=0.8 ,density=false)
    PyPlot.plt.xlabel("Distance CMD")
    PyPlot.plt.ylabel("N")
    PyPlot.plt.grid(true)
    PyPlot.plt.savefig(plotdir*"/"*plotfile)
    # PyPlot.plt.show()
end


## testing stars with very similar cmd
function(df::GaiaClustering.Df, dfcart::GaiaClustering.Df, dfnew::GaiaClustering.Df, dfcartnew::GaiaClustering.Df, idx)
    debug_red("entering testing tail...")
    doc= __transform_df(df, dfcart, idx)

    s=dfnew.ndata
    s1= length(dfnew.data[4,:])
    debug_red("$s $s1")
    dnew= __transform_df(dfnew, dfcartnew, 1:s)

    xcenter= median(doc.X) ; ycenter= median(doc.Y) ;  zcenter= median(doc.Z)
    dnew.X = dnew.X .- xcenter ; dnew.Y = dnew.Y .- ycenter ; dnew.Z = dnew.Z .- zcenter
    vlc= median(doc.vl) ; vbc= median(doc.vb) ; vrc= median(doc.vrad)

    radiusMax= 500  ## distance max to oc
    r2= radiusMax*radiusMax

    velocityMax= 10 ## velocity difference max
    v2= velocityMax*velocityMax

    dnewrad=filter(row -> (row.X*row.X+row.Y*row.Y+row.Y+row.Z*row.Z) < r2, dnew)
    dnewvel=filter(row -> ((row.vl .- vlc)*(row.vl .- vlc) + (row.vb .- vbc)*(row.vb .- vbc)) < v2, dnewrad)
    
    s= length(dfnew.data[1, :])
    srad= length(dnewrad.X)
    svel= length(dnewvel.X)

    debug_red("Filtering..")
    debug_red("$s $srad $svel")

    __plot_tail(dnew, doc, "test_tail1")
    __plot_tail(dnewrad, doc , "test_tail2")
    __plot_tail(dnewvel, doc , "test_tail3")

    doc= filter(:BmR0 => x -> !(ismissing(x) || isnothing(x) || isnan(x)), doc)
    dnewvel= filter(:BmR0 => x -> !(ismissing(x) || isnothing(x) || isnan(x)), dnewvel)

    idx , dist= __distance_cmd(doc,dnewvel)
    println(length(dnewvel.G))
    println(length(dist))

    ## cut with cmd
    cmdDistMax= 0.05
    idc= findall(x->(x< cmdDistMax),dist)
    __plot_dist_cmd(dist)
    s= length(idc)
    debug_red("CMD dist condition: $s")

    dnewcmd= dnewvel[idc,:]
    # println(idx[idc])
    __plot_tail(dnewcmd, doc,  "test_tail4")
    __density_count(dnewcmd.Y, dnewcmd.Z)
    
    
     a += 0
end

function __transform_df(df,dfcart, idx)
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
function __distance_cmd(df1, df2)
    n= length(df1.X)
    dfref= zeros(2,n)
    dfref[1,:] = df1.BmR0
    dfref[2,:] = df1.G

    kdtree = KDTree(dfref; leafsize = 20)

    npts= length(df2.G)
    pts= zeros(2,npts)
    pts[1,:] = df2.BmR0
    pts[2,:] = df2.G

    idxs, dists = nn(kdtree, pts)
    idx= collect(Iterators.flatten(idxs))
    dist= collect(Iterators.flatten(dists))
    idx=collect(1:npts)
    return(idx , dist)
end

function __density_count(xx, yy, nbin=100, xrange=[-100,100],yrange=[-100,100])
    data = (xx,yy) 
    stepx= (xrange[2]-xrange[1])/nbin
    stepy= (yrange[2]-yrange[1])/nbin 
    debug_red("density..")

    h = FHist.Hist2D((xx,yy), (xrange[1]:stepx:xrange[2], yrange[1]:stepy:yrange[2]))

    wav= atrous(h.hist.weights, 5)
    # wavFilt= thresholdingWav(wav,Normal())
    rec= addWav(wav,3,6)
    rec= permutedims(rec, [2, 1])
    nrec= size(rec)
    println(nrec)
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
    println(yv)
    println(yti)


    PyPlot.plt.figure(figsize=(9.0,8.0))
    ax= PyPlot.plt.subplot(1, 1, 1 )

    println("toto...")
    println(ax)
    ax.set_xticks(xv) ;  ax.set_xticklabels(xti)
    ax.set_yticks(yv) ;  ax.set_yticklabels(yti)

    vmin, vmax=  __level_dens(rec , 0.1, 5, norm="log")
    nlev= 20
    PyPlot.plt.contour(rec, nlev, vmin=vmin, vmax=vmax, linewidths= 0.2, colors= "black") 

    PyPlot.plt.savefig("test_density.png")

end

function __level_dens(dens,sigmin= 3, sigmax= 20, clip= 10)

    # sigma-clipping
    sigfirst= std(dens)
    mfirst= median(dens)
    sigfinal= std(dens[dens .< clip*sigfirst])
    vmin= sigmin*sigfinal
    vmax= sigmax*sigfinal

    return(vmin, vmax)
end

function __plot_surface_density(xx,yy, plotfile)

    r2d, dens, err2d= density2D(xx,yy,20)

    PyPlot.plt.figure(figsize=(9.0,8.0))
    ax= PyPlot.plt.subplot(1, 1, 1 )

    ax.set_yscale("log")
    ax.set_xscale("log")
    PyPlot.plt.xlim(r2d[1]*0.9, r2d[end]*1.1)
    PyPlot.plt.ylim(minimum(dens[dens .> 0])*0.5,maximum(dens)*1.5)
    PyPlot.plt.grid("on")
    PyPlot.plt.xlabel("radius (pc)")
    PyPlot.plt.ylabel("ρ")
    PyPlot.plt.scatter(r2d, dens , s=4, facecolor="blue" )
    PyPlot.plt.errorbar(r2d, dens, yerr=2 .* err2d, linewidth=0.5)


    PyPlot.plt.savefig(plotfile)
end