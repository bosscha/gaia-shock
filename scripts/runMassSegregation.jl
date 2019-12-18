## perform the estimation of the mass segregation parameters in the e2e_products directory
## create the plots for each cluster
## output the results in a .csv file

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
## Mass Segregation parameters#############################################
DTYPE= "2d"                 ## projection 2D/3D
DAVER= "geo"               ## average calculation: geo/ari
PERCMASSIVE= 20             ## percentile on the mag for massive/light stars

OUTFILE= "votlist-massSegregation.csv"
###########################################################################
############################functions #####################################
function plot_massSegregation(oc, figname, perc, kappa2d, kappa2derr)
    imass= select_massivestars(oc.gbar, perc, false)
    ilight= select_massivestars(oc.gbar, perc, true)

    xg= sum(oc.Y)/length(oc.Y) ;  yg= sum(oc.Z)/length(oc.Z)
    XX= oc.Y .- xg ; YY= oc.Z .- yg
    xmin= minimum(XX) ; xmax= maximum(XX)
    ymin= minimum(YY) ; ymax= maximum(YY)
    dl= max(xmax-xmin,ymax-ymin)*1.05 ; x2= (xmin+xmax)/2 ; y2= (ymin+ymax)/2

    fig= PyPlot.figure(figsize=(7,7))
    ax= PyPlot.subplot(111)
    PyPlot.grid("on")

    rr= sqrt.((XX .* XX) .+ (YY .* YY))
    siglight= std(rr[ilight]) ; sigmass= std(rr[imass])

    circle = PyPlot.plt.Circle((0,0), siglight, color="b", fill=false, linestyle="--")
    ax.add_artist(circle)
    circle = PyPlot.plt.Circle((0,0), sigmass, color="r", fill=false, linestyle="--")
    ax.add_artist(circle)

    PyPlot.title(figname[1:end-20])
    PyPlot.xlabel("Y(pc)") ; PyPlot.ylabel("Z(pc)")
    PyPlot.xlim([x2-dl/2,x2+dl/2]) ; PyPlot.ylim([y2-dl/2,y2+dl/2])

    PyPlot.plt.plot(XX[imass], YY[imass],"or",label="massive")
    PyPlot.plt.plot(XX[ilight], YY[ilight],"ob",label="light")
    PyPlot.legend()
    txt= @sprintf("kappa: %3.2f(%3.2f)",kappa2d,kappa2derr)
    PyPlot.text(x2-dl*0.45,y2-dl*0.45,txt)
    PyPlot.savefig(figname)
end
## main loop
    cd(ocdir)
    files= glob("*csv")
    cd(wdir)

    ilabel= []
    Κms2d= Array{Float64,1}(undef,0) ; Κms3d = Array{Float64,1}(undef,0) ; k3=  Array{Float64,1}(undef,0)
    errΚms2d= Array{Float64,1}(undef,0) ; errΚms3d= Array{Float64,1}(undef,0) ; errk3= Array{Float64,1}(undef,0)
    votname= Array{String,1}(undef,0)

    i= 0
    for f in files
        i+= 1
        oc= CSV.read("$ocdir/$f" , delim= ";")

        vot= @sprintf("%s.vot", f[1:end-7])
        println(vot)
        push!(votname, vot )

        y, yerr= get_kappaMS(oc[:gbar],oc[:X], oc[:Y] , oc[:Z],PERCMASSIVE , "3d", DAVER, false)
        push!(ilabel, i)
        push!(Κms3d, y)
        push!(errΚms3d, yerr)

        y, yerr= get_kappaMS(oc[:gbar],oc[:X], oc[:Y] , oc[:Z], PERCMASSIVE, "2d", DAVER, false)
        push!(Κms2d, y)
        push!(errΚms2d, yerr)

        y, yerr= get_kappaMS(oc[:gbar],oc[:vl], oc[:vl] , oc[:vb], PERCMASSIVE, "2d", DAVER, false)
        push!(k3, y)
        push!(errk3, yerr)

        cd(plotdir)
        figname= @sprintf("%s-massSegregation.png", f[1:end-7])
        plot_massSegregation(oc, figname, PERCMASSIVE, Κms2d[end],errΚms2d[end])
        cd(wdir)
    end

    df= DataFrames.DataFrame()
    df.votname= votname
    df.kappa2d= Κms2d ; df.kappa2derr= errΚms2d
    df.kappa3d= Κms3d ; df.kappa3derr= errΚms3d
    df.kappavel= k3 ; df.kappavelerr= errk3

    CSV.write(OUTFILE, df, delim=";")
    println(df)

println(median(Κms2d))
println(std(Κms2d))
println(median(Κms3d))
println(std(Κms3d))
println(median(k3))
println(std(k3))


end
