## estimae the stellar mass using a known age of the ]OC
##
rootdir = ENV["GAIA_ROOT"]
push!(LOAD_PATH,"$rootdir/master/src")
using GaiaClustering

using PyCall
## PopStar modules...
synthetic= pyimport("popstar.synthetic")
evolution= pyimport("popstar.evolution")
atmospheres= pyimport("popstar.atmospheres")
reddening= pyimport("popstar.reddening")

using Distributions , Statistics , DataFrames , CSV
using Distances, Random , Printf, Glob

import PyPlot

## directory
wdir=    "$rootdir/e2e_products"
plotdir= "$wdir/plotsIndividuals"
ocdir=   "$wdir/oc"
isodir=  "$wdir/isochrones"

let
## ISOCHRONES parameters
VOTAGE= "votlist-age.csv"    ## csv file with votname,log_age

## Main Loop ####################################################
cd(ocdir)
files= glob("*-oc.csv")
cd(wdir)
votage= CSV.read(VOTAGE, delim=";")

for f in files
    println("## $f")
    fvot= f[1:end-7] ; vname= "$fvot.vot" ; fname= f[1:end-4]
    ocmassfile= "$ocdir/$fname-mass.csv"
    logAge= votage[(votage.votname .== vname), :log_age][1]
    println("$ocdir/$f")
    df= CSV.read("$ocdir/$f" , delim= ";")

    dist= median(df.distance)
    AKs= 0.0                        # extinction in mags
    metallicity= 0.0                # Metallicity in [M/H]
    nstar= length(df.distance)

    ## Define evolution/atmosphere models and extinction law
    if logAge>6.5
        evo_model = evolution.MergedBaraffePisaEkstromParsec()
    else
        evo_model= evolution.MISTv1()
    end

    atm_func = atmospheres.get_merged_atmosphere
    red_law = reddening.RedLawHosek18b()

    ## Filter
    filt_list = [ "gaia,dr2_rev,Gbp", "gaia,dr2_rev,G"  ,"gaia,dr2_rev,Grp"]

    #####################################
    println("### Computing isochrones..")
    iso = synthetic.IsochronePhot(logAge, AKs, dist, metallicity=metallicity,
                            evo_model=evo_model, atm_func=atm_func,
                            red_law=red_law, filters=filt_list, iso_dir= isodir)

    data= iso.points
    ndata= size(data)[1]

    mass= zeros(ndata)
    Gmag= zeros(ndata)
    BPmag= zeros(ndata)
    RPmag= zeros(ndata)
    color= zeros(ndata)

    for i in 1:ndata
        mass[i]= data[i][4]
        Gmag[i]= data[i][10]
        BPmag[i]= data[i][9]
        RPmag[i]= data[i][11]
        color[i]= BPmag[i]-RPmag[i]
    end

    ## Cluster
    BmR= df.bp .- df.rp
    GMAG= df.gbar + 5 .* log10.(df.distance) .- 17.

    ### Stellar mass
    ################
    nxy= length(mass)
    A= Array{Float64}(undef,2,nxy)
    for i in 1:nxy
        A[1,i]= color[i]
        A[2,i]= Gmag[i]
    end

    d= Distances.Euclidean()
    mstar= zeros(length(BmR))

    for i in 1:length(BmR)
        starm= [BmR[i]; GMAG[i]]
        r= Distances.colwise(d, A, starm)

        idx= argmin(r)
        mstar[i]= mass[idx]
    end

    df.mass= mstar

    ### save the OC df
    CSV.write(ocmassfile, df, delim=";")
    println("\n\n")
end
###
end
