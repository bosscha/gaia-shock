### Functions for plotting stellar cluster results
###

#### aux. function to plot text in a subplot
## plot array of text in a box
function show_text(posx,pys, text, ywidth = 1.)
    posy = pys
    dy = ywidth / length(text)
    for t in text
        PyPlot.plt[:text](posx,posy,t)
        posy += dy
    end
end

### plot the stats of the DBSCAN parameters from the MCMC
function plot_dbscan_mcmc(plotdir, voname, mc::GaiaClustering.mc , showplot = true)
    PyPlot.plt[:figure](figsize=(10.0,12.0))

    nbins = 50
    PyPlot.plt[:subplot](3, 2, 1  )
    h = PyPlot.plt[:hist](mc.eps,nbins, alpha=0.75)
    PyPlot.plt[:xlabel]("ϵ")
    PyPlot.plt[:grid](true)
    
    nbins = 20
    PyPlot.plt[:subplot](3, 2, 2 )
    h = PyPlot.plt[:hist](mc.mne,nbins, alpha=0.75)
    PyPlot.plt[:xlabel]("min_neighbor")  
    PyPlot.plt[:grid](true)
    
    nbins = 20
    PyPlot.plt[:subplot](3, 2, 3 )
    h = PyPlot.plt[:hist](mc.mcl,nbins, alpha=0.75)
    PyPlot.plt[:xlabel]("min_cluster") 
    PyPlot.plt[:grid](true)
    
    ### text 

    PyPlot.plt[:subplot](3, 2, 4)
    PyPlot.plt[:axis]("off")
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
    PyPlot.plt[:subplot](3, 2, 6 )
    PyPlot.plt[:axis]("on")
    PyPlot.plt[:axis]("on")
    h = PyPlot.plt[:hist](mc.qn,nbins, color = "g", alpha=0.75)
    PyPlot.plt[:xlabel]("Qn")     
    PyPlot.plt[:grid](true)
    
    nbins = 50
    PyPlot.plt[:subplot](3, 2, 5 )
    h = PyPlot.plt[:hist](mc.qc,nbins, color="g", alpha=0.75)
    PyPlot.plt[:grid](true)

    PyPlot.plt[:xlabel]("Qc") 
    figname = plotdir*"/"*voname*"-dbscanMCMC.png"
    PyPlot.plt[:savefig](figname)
    if showplot PyPlot.plt[:show]() end
end
############
