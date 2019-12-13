## perform the estimation of the mass segregation parameters in the e2e_products directory
## create the plots for each cluster
## output the results in a .csv file

rootdir= ENV["GAIA_ROOT"]
push!(LOAD_PATH,"$rootdir/master/src")
using GaiaClustering
using Printf, Glob, Statistics

import PyPlot , CSV

## directory
wdir    = "$rootdir/e2e_products"
plotdir = "$wdir/plotIndividuals"
ocdir=    "$wdir/oc"

cd(wdir)

let
## Mass Segregation parameters#############################################
DTYPE= "2d"                 ## projection 2D/3D
DAVER= "geo"               ## average calculation: geo/ari
PERCMASSIVE= 15             ## percentile on the mag for massive/light stars

###########################################################################
## main loop
    cd(ocdir)
    files= glob("*csv")
    cd(wdir)

    ilabel= []
    Κms2d= Array{Float64,1}(undef,0) ; Κms3d = Array{Float64,1}(undef,0) ; k3=  Array{Float64,1}(undef,0)
    errΚms2d= Array{Float64,1}(undef,0) ; errΚms3d= Array{Float64,1}(undef,0) ; errk3= Array{Float64,1}(undef,0)

    i= 0
    for f in files
        i+= 1
        oc= CSV.read("$ocdir/$f" , delim= ";")

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
    end

println(median(Κms2d))
println(std(Κms2d))
println(median(Κms3d))
println(std(Κms3d))
println(median(k3))
println(std(k3))
end


############################functions ##########################################
function plot_massSegregation()
end
