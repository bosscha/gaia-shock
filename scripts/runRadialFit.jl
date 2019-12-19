## perform the Cauchy radial fit on the ocs
##

rootdir= ENV["GAIA_ROOT"]
push!(LOAD_PATH,"$rootdir/master/src")
using GaiaClustering
using Printf, Glob, Statistics

import PyPlot , CSV , DataFrames

## directory
wdir    = "$rootdir/e2e_products"
plotdir = "$wdir/plotsIndividuals"
ocdir=    "$wdir/oc"

cd(wdir)

let

## fit parameters##########################################################
       NTEST= 100 ; NITER= 30000 ; NBURNIN= 5000 ##MCMC parameters for the fit
       OUTFILE= "votlist-CauchyFit.csv"
###########################################################################
    @printf("## Compute the Cauchy radial fit on the OCs\n")
    @printf("##\n")
############################functions #####################################
function plotFitDensityFig(filelist, θlist, ifig, nrow, ncol, ocdir)
    fig= PyPlot.figure(figsize=(15,10))

    for i in 1:length(filelist)
        ocname= filelist[i][1:end-4]
        θ= θlist[i]
        ocfile= "$ocdir/$ocname-oc.csv"
        plotDensity(ocfile, θ,nrow, ncol, i, ocname)
    end

    figname= "densityOC-$ifig.png"
    PyPlot.savefig(figname)
    PyPlot.close(fig)
end


function plotDensity(ocfile, θ, irow, icol, nplot, title)
    oc= CSV.read(ocfile, delim= ";")
    # binned density
    nstar= length(oc.X)
    nbin= min(trunc(Int, nstar/10),20)
    r2d,ρ2d,err2d= density2D(oc.Y, oc.Z, nbin)
    ρ2dfit= model_rad(r2d, θ, fdens1)

    ax= PyPlot.subplot(irow, icol, nplot)
    ax.set_xlabel("r [pc]") ; ax.set_ylabel("ρ")
    ax.set_yscale("log") ; ax.set_xscale("log")
    ax.set_xlim(r2d[1]*0.9, r2d[end]*1.1)
    ax.set_ylim(minimum(ρ2d[ρ2d .> 0])*0.3,maximum(ρ2d)*1.5)
    PyPlot.grid("on")
    PyPlot.scatter(r2d, ρ2d , s=4, facecolor="blue" )
    PyPlot.errorbar(r2d, ρ2d, yerr=2 .* err2d, linewidth=0.5)
    PyPlot.plot(r2d, ρ2dfit, "k-", linewidth=1)

    ## TEXT
    ax.text(0.02, 0.01, title, verticalalignment= "bottom" , horizontalalignment="left",
        transform=ax.transAxes, color="green", fontsize=15)
    val=θ.C ; txt= @sprintf("C= %3.3f",val)
    ax.text(0.02, 0.2, txt, verticalalignment= "bottom" , horizontalalignment="left",
        transform=ax.transAxes, color="red", fontsize=12)
    val=θ.s ; txt= @sprintf("s= %3.3f [pc]",val)
    ax.text(0.02, 0.3, txt, verticalalignment= "bottom" , horizontalalignment="left",
        transform=ax.transAxes, color="red", fontsize=12)
    val=θ.m ; txt= @sprintf("m= %3.3f",val)
    ax.text(0.02, 0.4, txt, verticalalignment= "bottom" , horizontalalignment="left",
        transform=ax.transAxes, color="red", fontsize=12)
    val= nstar ; txt= @sprintf("Nstar= %d",val)
    ax.text(0.02, 0.52, txt, verticalalignment= "bottom" , horizontalalignment="left",
        transform=ax.transAxes, color="blue", fontsize=12)
end

## main loop
    cd(ocdir)
    files= glob("*csv")
    cd(wdir)

    df = DataFrames.DataFrame(votname= String[], C= Float64[], Cerr= Float64[], s= Float64[], serr= Float64[], m= Float64[],
merr= Float64[])
    i= 0 ; noc= 0 ; nfail= 0
    for f in files
        i+= 1
        ocfile= "$ocdir/$f"
        oc= CSV.read("$ocdir/$f" , delim= ";")
        println("### ocfile: $f")
        println("### fitting the radial profile....")
        votname= @sprintf("%s.vot", f[1:end-7])

        nstar= length(oc.X)
        nbin= min(trunc(Int, nstar/10),20)

        verbose= false
        θfit, θfiterr, fitfound= spatialParameter(ocfile, NTEST, nbin, NITER, verbose, NBURNIN)
        println("### Fit:")
        println(θfit)
        println(θfiterr)

        push!(df, [votname,θfit.C , θfiterr.C, θfit.s , θfiterr.s, θfit.m , θfiterr.m])

        noc += 1
        if !fitfound
            println("################################ Failing...########### ")
            println("###################################################### ")
            nfail += 1
        end
        println("## OCs analyzed: $noc")
        println("## Fit failed: $nfail \n\n")
    end
    CSV.write(OUTFILE, df, delim=";")

    ######################## plots
    println("## Plotting the individual fits...")
    cd(plotdir)
    rowperfig= 4 ; colperfig= 3 ; plotperfig= rowperfig*colperfig

    ifig= 1 ; plotaccum=0
    filevot= [] ; θlist= [] ; θerrlist= []
    for item in 1:length(df.votname)
        push!(filevot, df.votname[item])
        θ= GaiaClustering.modelCauchy(df.C[item], df.s[item], df.m[item])
        push!(θlist, θ)
        θerr= GaiaClustering.modelCauchy(df.Cerr[item], df.serr[item], df.merr[item])
        push!(θerrlist, θerr)
        plotaccum += 1

        if plotaccum == plotperfig
            plotFitDensityFig(filevot, θlist, ifig,rowperfig, colperfig , ocdir)
            println("## Fig $ifig done.")
            filevot= [] ; θlist= [] ; θerrlist= []
            plotaccum= 0 ; ifig += 1
        end
    end
    if plotaccum > 0
        plotFitDensityFig(filevot, θlist, ifig,rowperfig, colperfig , ocdir)
        println("## Fig $ifig done.")
    end
    cd(wdir)
end
