### Functions for plotting stellar cluster results
###

#### aux. function to plot text in a subplot
## plot array of text in a box
function show_text(posx,posy, text, ywidth = 1.)
    dct= Dict("color"=> "black", "fontsize"=> 11)
    dy = ywidth / length(text)
    for t in text
        PyPlot.plt.text(posx,posy,t, dct)
        posy += dy
    end
end

### plot the stats of the DBSCAN parameters from the MCMC
function plot_dbscan_mcmc(plotdir, voname, mc::GaiaClustering.mc , showplot = true)
    PyPlot.plt[:figure](figsize=(10.0,12.0))

    nbins = 50
    PyPlot.plt.subplot(3, 2, 1  )
    h = PyPlot.plt.hist(mc.eps,nbins, alpha=0.75)
    PyPlot.plt.xlabel("ϵ")
    PyPlot.plt.grid(true)

    nbins = 20
    PyPlot.plt.subplot(3, 2, 2 )
    h = PyPlot.plt.hist(mc.mne,nbins, alpha=0.75)
    PyPlot.plt.xlabel("min_neighbor")
    PyPlot.plt.grid(true)

    nbins = 20
    PyPlot.plt.subplot(3, 2, 3 )
    h = PyPlot.plt.hist(mc.mcl,nbins, alpha=0.75)
    PyPlot.plt.xlabel("min_cluster")
    PyPlot.plt.grid(true)

    ### text

    PyPlot.plt.subplot(3, 2, 4)
    PyPlot.plt.axis("off")
    ## text to display
        text =[]
        push!(text, "VOT: "*voname)
        vtext  = mean(mc.eps)
        v2text = std(mc.eps)
        txt = "ϵ : $vtext +/- $v2text "
        push!(text,txt)
        vtext  = mean(mc.mne)
        v2text = std(mc.mne)
        txt = "min_neigh : $vtext +/- $v2text "
        push!(text,txt)
        vtext  = mean(mc.mcl)
        v2text = std(mc.mcl)
        txt = "min_clus : $vtext +/- $v2text "
        push!(text,txt)
        vtext  = mean(mc.qn)
        v2text = std(mc.qn)
        txt = "Qn : $vtext +/- $v2text "
        push!(text,txt)
        vtext  = mean(mc.qc)
        v2text = std(mc.qc)
        txt = "Qc : $vtext +/- $v2text "
        push!(text,txt)
        show_text(-0.01,0.0, text)

    #######
    nbins = 50
    PyPlot.plt.subplot(3, 2, 6 )
    PyPlot.plt.axis("on")
    PyPlot.plt.axis("on")
    h = PyPlot.plt.hist(mc.qn,nbins, color = "g", alpha=0.75)
    PyPlot.plt.xlabel("Qn")
    PyPlot.plt.grid(true)

    nbins = 50
    PyPlot.plt.subplot(3, 2, 5 )
    h = PyPlot.plt.hist(mc.qc,nbins, color="g", alpha=0.75)
    PyPlot.plt.grid(true)

    PyPlot.plt.xlabel("Qc")
    figname = plotdir*"/"*voname*".mcmc.png"
    PyPlot.plt.savefig(figname)
    if showplot PyPlot.plt.show() end
end


### plot the stats of the DBSCAN parameters from the MCMC
function plot_dbscanfull_mcmc(plotdir, voname, mc::mcfull , showplot = true)
    patch= pyimport("matplotlib.patches")
    PyPlot.plt.figure(figsize=(12.0,13.0))

    nbins = 50
    PyPlot.plt.subplot(3, 3, 1  )
    h = PyPlot.plt.hist(mc.eps,nbins, alpha=0.75)
    PyPlot.plt.xlabel("ϵ")
    PyPlot.plt.grid(true)

    nbins = 20
    PyPlot.plt.subplot(3, 3, 2 )
    h = PyPlot.plt.hist(mc.mne,nbins, alpha=0.75)
    PyPlot.plt.xlabel("min_neighbor")
    PyPlot.plt.grid(true)

    nbins = 20
    PyPlot.plt.subplot(3, 3, 3 )
    h = PyPlot.plt.hist(mc.mcl,nbins, alpha=0.75)
    PyPlot.plt.xlabel("min_cluster")
    PyPlot.plt.grid(true)

    nbins = 50
    PyPlot.plt.subplot(3, 3, 4 )
    h = PyPlot.plt.hist(mc.w3d,nbins, alpha=0.75)
    PyPlot.plt.xlabel("W3d")
    PyPlot.plt.grid(true)

    nbins = 50
    PyPlot.plt.subplot(3, 3, 5 )
    h = PyPlot.plt.hist(mc.wvel,nbins, alpha=0.75)
    PyPlot.plt.xlabel("Wvel")
    PyPlot.plt.grid(true)

    nbins = 50
    PyPlot.plt.subplot(3, 3, 6 )
    h = PyPlot.plt.hist(mc.whrd,nbins, alpha=0.75)
    PyPlot.plt.xlabel("Whrd")
    PyPlot.plt.grid(true)


    #######
    nbins = 50
    PyPlot.plt.subplot(3, 3, 7 )
    PyPlot.plt.axis("on")
    PyPlot.plt.axis("on")
    h = PyPlot.plt.hist(mc.qn,nbins, color = "g", alpha=0.75)
    PyPlot.plt.xlabel("Qn")
    PyPlot.plt.grid(true)

    nbins = 50
    PyPlot.plt.subplot(3, 3, 8 )
    h = PyPlot.plt.hist(mc.qc,nbins, color="g", alpha=0.75)
    PyPlot.plt.xlabel("Qc")
    PyPlot.plt.grid(true)

        ### text

    axt= PyPlot.plt.subplot(3, 3, 9)
    PyPlot.plt.axis("off")
    ## text to display
        dct= Dict( :color=> "black", :fontsize=> 10)
        text =[]
        push!(text, "votable: "*voname)
        vtext  =  fmt("3.3f", median(mc.eps))
        v2text =  fmt("3.3f", std(mc.eps))
        txt = "ϵ : $vtext +/- $v2text "
        push!(text,txt)
        vtext  =  fmt("3.3f", median(mc.mne))
        v2text =  fmt("3.3f", std(mc.mne))
        txt = "min_neigh : $vtext +/- $v2text "
        push!(text,txt)
        vtext  =  fmt("3.3f", median(mc.mcl))
        v2text =  fmt("3.3f", std(mc.mcl))
        txt = "min_clus : $vtext +/- $v2text "
        push!(text,txt)
        vtext  =  fmt("3.3f", median(mc.w3d))
        v2text =  fmt("3.3f", std(mc.w3d))
        txt = "W3d : $vtext +/- $v2text "
        push!(text,txt)
        vtext  =  fmt("3.3f", median(mc.wvel))
        v2text =  fmt("3.3f", std(mc.wvel))
        txt = "Wvel : $vtext +/- $v2text "
        push!(text,txt)
        vtext  =  fmt("3.3f", median(mc.whrd))
        v2text =  fmt("3.3f", std(mc.whrd))
        txt = "Whrd : $vtext +/- $v2text "
        push!(text,txt)
        vtext  =  fmt("3.3f", median(mc.qn))
        v2text =  fmt("3.3f", std(mc.qn))
        txt = "Qn : $vtext +/- $v2text "
        push!(text,txt)
        vtext  =  fmt("3.3f", median(mc.qc))
        v2text =  fmt("3.3f", std(mc.qc))
        txt = "Qc : $vtext +/- $v2text "
        push!(text,txt)
        show_text(-0.01,0.0, text)

        rec= patch.Rectangle((-0.07, -0.05), 1.2 , 1.05, color="khaki",
            alpha= 0.4, clip_on=false)
        axt.add_artist(rec)

    PyPlot.plt.xlabel("Qc")
    figname = plotdir*"/"*voname*".mcmc.png"
    PyPlot.plt.savefig(figname)
    if showplot PyPlot.plt.show() end
end


function plot_cluster(plotdir, voname, indx, sc::GaiaClustering.SCproperties, df::GaiaClustering.Df, showplot = true , cmap = "gist_stern")

    println("### Cluster plot is centered in Y,Z...")

    PyPlot.plt.figure(figsize=(13.0,12.0))
    PyPlot.plt.subplot(3, 3, 1 , xlim = [-20,20] , ylim = [-20,20])

    xx = df.data[2,indx] .- mean(df.data[2,indx])
    yy = df.data[3,indx] .- mean(df.data[3,indx])

    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("Y (pc)")
    PyPlot.plt.ylabel("Z (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(3, 3, 2 , ylim = [-20,20])
    xx = df.data[1,indx]
    yy = df.data[3,indx] .- mean(df.data[3,indx])
    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("X (pc)")
    PyPlot.plt.ylabel("Z (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(3, 3, 4 , xlim = [-20,20])
    xx = df.data[2,indx] .- mean(df.data[2,indx])
    yy = df.data[1,indx]
    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("Y (pc)")
    PyPlot.plt.ylabel("X (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(3, 3, 3 )
    xx = df.data[1,indx]
    yy = df.raw[13,indx]
    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("X(pc)")
    PyPlot.plt.ylabel("Vrad (km/s)")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(3, 3, 5)
    PyPlot.plt.axis("off")
    ## text to display
    text =[]
    v = sc.nstars ; txt = "N stars   : $v" ; push!(text,txt)
    v = fmt("3.1f",sc.distance) ; txt = "Distance  : $v (pc)" ; push!(text,txt)
    v = fmt("3.3f",sc.l) ; txt = "l         : $v (degree)" ; push!(text,txt)
    v = fmt("3.3f",sc.b) ; txt = "b         : $v (degree)" ; push!(text,txt)
    v = fmt("3.2f", sc.vl) ; txt = "Vl      : $v (km/s)" ; push!(text,txt)
    v = fmt("3.2f",sc.vb) ; txt = "Vb      : $v (km/s)"; push!(text,txt)
    v = fmt("3.2f",sc.vrad) ; txt = "Vradial   : $v (km/s)"; push!(text,txt)
    v = fmt("3.2f",sc.xdisp) ; txt = "X disp.   : $v (pc)" ; push!(text,txt)
    v = fmt("3.2f",sc.ydisp) ; txt = "Y disp.   : $v (pc)" ; push!(text,txt)
    v = fmt("3.2f",sc.zdisp) ; txt = "Z disp.   : $v (pc)" ; push!(text,txt)
    v = fmt("3.2f",sc.vldisp) ; txt = "Vl disp. : $v (km/s)" ; push!(text,txt)
    v = fmt("3.2f",sc.vbdisp) ; txt = "Vb disp.: $v (km/s)" ; push!(text,txt)
    v = fmt("3.2f",sc.vraddisp) ; txt = "Vradial disp.: $v (km/s)" ; push!(text,txt)
    show_text(-0.01,0.0, text , 1.0)

    PyPlot.plt.subplot(3, 3, 7 )
    PyPlot.plt.axis("on")
    xx = df.data[7,indx]
    yy = -df.data[6,indx]
    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("G-Rp")
    PyPlot.plt.ylabel("G")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(3, 3, 8 )
    xx = df.data[4,indx]
    yy = df.data[5,indx]
    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("Vl (km/s)")
    PyPlot.plt.ylabel("Vb (km/s)")
    PyPlot.plt.grid(true)


    figname = plotdir*"/"*voname*".cluster.png"
    PyPlot.plt.savefig(figname)
    if showplot PyPlot.plt.show() end
end

function plot_cluster2(plotdir, voname, indx, sc::GaiaClustering.SCproperties2, df::GaiaClustering.Df,
    showplot = true , extra= 0, cmap = "gist_stern")
    patch= pyimport("matplotlib.patches")

    println("### Cluster plot is centered in Y,Z...")
    PyPlot.plt.rcParams["font.size"]= 25

    PyPlot.plt.figure(figsize=(13.0,12.0))
    PyPlot.plt.subplot(3, 3, 1 , xlim = [-20,20] , ylim = [-20,20])

    xx = df.data[2,indx] .- mean(df.data[2,indx])
    yy = df.data[3,indx] .- mean(df.data[3,indx])

    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("Y (pc)")
    PyPlot.plt.ylabel("Z (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(3, 3, 2 , ylim = [-20,20])
    xx = df.data[1,indx]
    yy = df.data[3,indx] .- mean(df.data[3,indx])
    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("X (pc)")
    PyPlot.plt.ylabel("Z (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(3, 3, 4 , xlim = [-20,20])
    xx = df.data[2,indx] .- mean(df.data[2,indx])
    yy = df.data[1,indx]
    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("Y (pc)")
    PyPlot.plt.ylabel("X (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(3, 3, 3 )
    xx = df.data[1,indx]
    yy = df.raw[13,indx]
    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("X(pc)")
    PyPlot.plt.ylabel("Vrad (km/s)")
    PyPlot.plt.grid(true)

    axt= PyPlot.plt.subplot(3, 3, 5)
    PyPlot.plt.axis("off")

    ## text to display
    dct= Dict("color"=> "black", "fontsize"=> 10)
    text =[]
    v = sc.nstars ; txt = "N stars : $v" ; push!(text,txt)
    v = fmt("3.1f",sc.distance) ; txt = "Distance : $v (pc)" ; push!(text,txt)
    v1 = fmt("3.3f",sc.l) ; v2 = fmt("3.3f",sc.b) ;
    txt = "l , b : $v1  ,  $v2  (degree)" ; push!(text,txt)
    v1 = fmt("3.3f",sc.ra) ; v2 = fmt("3.3f",sc.dec) ;
    txt = "RA , Dec : $v1  ,  $v2  (degree)" ; push!(text,txt)

    v1 = fmt("3.3f",sc.vl) ; v2 = fmt("3.3f",sc.vb) ;
    txt = "vl , vb : $v1  ,  $v2  (km/s)" ; push!(text,txt)
    v = fmt("3.2f",sc.vrad) ; txt  = "Vradial : $v (km/s)"; push!(text,txt)
    v = fmt("3.2f",sc.xdisp) ; txt = "X disp. : $v (pc)" ; push!(text,txt)
    v = fmt("3.2f",sc.ydisp) ; txt = "Y disp. : $v (pc)" ; push!(text,txt)
    v = fmt("3.2f",sc.zdisp) ; txt = "Z disp. : $v (pc)" ; push!(text,txt)
    v = fmt("3.2f",sc.vldisp) ; txt = "Vl disp. : $v (km/s)" ; push!(text,txt)
    v = fmt("3.2f",sc.vbdisp) ; txt = "Vb disp. : $v (km/s)" ; push!(text,txt)
    v = fmt("3.2f",sc.vraddisp) ; txt = "Vradial disp. : $v (km/s)" ; push!(text,txt)
    show_text(-0.01,-0.1, text , 1.1 )

    if extra != 0
        text =[]
        v1= "$(extra.votname[1])"
        txt = "Votable : $v1" ; push!(text,txt)
        v1= "$(extra.cycle[1])"
        txt = "Cycle : $v1" ; push!(text,txt)
        v1= fmt("3.3f",extra.qc[1])
        txt = "Qc : $v1" ; push!(text,txt)
        v1= fmt("3.3f",extra.score_cycle[1])
        txt = "Score : $v1" ; push!(text,txt)
        v1= fmt("3.3f",sc.offdeg)
        txt = "Offset : $v1 (degree)" ; push!(text,txt)
        v1= fmt("3.3f",sc.edgratg)
        txt = "Edge ratio(g) : $v1 " ; push!(text,txt)
        v1= fmt("3.3f",sc.edgratm)
        txt = "Edge ratio(m) : $v1 " ; push!(text,txt)
        show_text(1.2,-0.1, text , 0.64 )

        rec= patch.Rectangle((-0.07, -0.15), 2.3, 1.15, color="salmon", alpha= 0.4, clip_on=false)
        axt.add_artist(rec)
    end

    PyPlot.plt.subplot(3, 3, 7 )
    PyPlot.plt.axis("on")
    xx = df.data[7,indx]
    yy = -df.data[6,indx]
    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("G-Rp")
    PyPlot.plt.ylabel("G")
    PyPlot.plt.grid(true)

    PyPlot.plt.subplot(3, 3, 8 )
    xx = df.data[4,indx]
    yy = df.data[5,indx]
    PyPlot.plt.scatter(xx, yy , s = 1.0 )
    PyPlot.plt.xlabel("Vl (km/s)")
    PyPlot.plt.ylabel("Vb (km/s)")
    PyPlot.plt.grid(true)


    figname = plotdir*"/"*voname*".cluster.png"
    PyPlot.plt.savefig(figname)
    if showplot PyPlot.plt.show() end
end
