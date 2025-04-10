## functions to perform isochrones fitting on the Gaia CMD 
##
## load the isochrone models
## Read a MIST cmd file
## Read a MIST cmd file
function mist_df(file)
    start=false
    ntrack= 0 ; df=0
    arrdf= []

    for line in eachline(file)
         if !isempty(line)
            line=lstrip(line)
            s1= line[1]

            if start && s1 != '#'
                sl= split(line)

                log10_isochrone_age_yr= tryparse(Float64,sl[2])
                initial_mass= tryparse(Float64,sl[3])
                star_mass=  tryparse(Float64,sl[4])
                feh_init = tryparse(Float64,sl[8])
                feh= tryparse(Float64,sl[9])
                Gaia_G_EDR3= tryparse(Float64,sl[31])
                Gaia_BP_EDR3=  tryparse(Float64,sl[32])
                Gaia_RP_EDR3= tryparse(Float64,sl[33])

                Gaia_BPmRP_EDR3= Gaia_BP_EDR3 - Gaia_RP_EDR3

                push!(df,(log10_isochrone_age_yr, initial_mass, star_mass, feh_init, feh, Gaia_G_EDR3,Gaia_BP_EDR3,Gaia_RP_EDR3,Gaia_BPmRP_EDR3))

            elseif start && s1 == '#'
                    start= false
                    ## new track ......
                    ntrack += 1
                    if ntrack > 0
                        push!(arrdf,df)
                    end
            end

            if !start && s1 != '#'
                start= true
                sl= split(line)

                log10_isochrone_age_yr= tryparse(Float64,sl[2])
                initial_mass= tryparse(Float64,sl[3])
                star_mass=  tryparse(Float64,sl[4])
                feh_init = tryparse(Float64,sl[8])
                feh= tryparse(Float64,sl[9])
                Gaia_G_EDR3= tryparse(Float64,sl[31])
                Gaia_BP_EDR3=  tryparse(Float64,sl[32])
                Gaia_RP_EDR3= tryparse(Float64,sl[33])

                Gaia_BPmRP_EDR3= Gaia_BP_EDR3 - Gaia_RP_EDR3

                df=DataFrame(log10_isochrone_age_yr=log10_isochrone_age_yr, initial_mass=initial_mass, star_mass=star_mass, feh_init=feh_init, feh=feh,
                             Gaia_G_EDR3=Gaia_G_EDR3,Gaia_BP_EDR3=Gaia_BP_EDR3,Gaia_RP_EDR3=Gaia_RP_EDR3, Gaia_BPmRP_EDR3=Gaia_BPmRP_EDR3)

            end
         end
    end
    println("### $ntrack MIST isochrones read ...")
    return(arrdf)
end

################################################################
#Read a collection of MIST files and concat them
function read_isochrones(f)
    arrIso=[]
    for file in f
        debug_red("# Reading MIST: $file")
        arr= mist_df(file)
        append!(arrIso,arr)
    end
    return(arrIso)
end
####################################################################
### Updating magnitude of a solution for an oc (csv format)
## filt: remove row with undefined values
function update_mag(df,filt=true)
    df.gbar0 = df.gbar .- df.ag
    df.BmR0 = df.bp .- df.rp .- df.ebmr
    ## absolute magnitude
    df.G = df.gbar0 .+ 5 .- 5log10.(df.distance)

    debug_red("G magnitude corrected for extinction...")
    if filt
        df= filter(:G => G -> !any(f -> f(G), (ismissing, isnothing, isnan)), df)
    end
    return(df)
end
#####################################################################
## weight for fitting isochrones with more wgt for larger magnitude (less pts)
## df: dataframe for the oc with the updated magnitude (update_mag)
## step: bin in magnitude for Histogram
## minpts: minimumpts in the bin to consider for weight
function  weight_cmd(df, step=0.25, minpts=1)
    debug_red("Weight...")

    magmin=minimum(df.G)
    magmax= maximum(df.G)
    # step=(magmax-magmin)/nbins
    h=  fit(Histogram, df.G , magmin:step:magmax)

    edg= collect(h.edges[1])
    wgt=zeros(length(df.G))

    nedg= length(edg)-1
    for i in 1:nedg
        ind= findall(x-> edg[i]<=x<edg[i+1], df.G)

        # minimum pts to consider..
        if h.weights[i] > minpts
            wgt[ind]  .= 1/ h.weights[i]
        end
    end

    sumwgt= sum(wgt)
    wgt = wgt ./ sumwgt
    return(wgt)
end
#############################################
## total distance of pts in CMD to isochrone
## df: oc solution
## iso: isochrone df (MIST)
## wgt: weight for each pt
function dist_cmd2iso(df, iso, wgt)
    n= length(iso.Gaia_G_EDR3)

    diso= zeros(2,n)
    diso[1,:] = iso.Gaia_BPmRP_EDR3
    diso[2,:] = iso.Gaia_G_EDR3

    # debug_red("dist_cmd ...........")

    kdtree = KDTree(diso; leafsize = 10)

    npts= length(df.G)
    pts= zeros(2,npts)
    pts[1,:] = df.BmR0
    pts[2,:] = df.G

    idxs, dists = knn(kdtree, pts, 1, true)

    dists = dists .* wgt
    total= sum(dists)

    # debug_red("dists total ... $total")

    return(total)
end
################################################################################
#### main function to fit the best isochrone
## df: OC solution
## arrIso: array of the isochrones to be tested
## wgt: weight for the pts in the OC solution
function fit_isochrone(df,arrIso,wgt)
    println("## starting isochrone fitting ...")

    println("### removing rows with NaN for fitting... ")
    df= filter(:G => G -> !any(f -> f(G), (ismissing, isnothing, isnan)), df)

    idiso= 1 ; distmin= 1e9
    bestiso= [] ; bestage= 0; bestfeh= 0

    for iso in arrIso
        age= iso.log10_isochrone_age_yr[1]
        feh= iso.feh[1]

        dist= dist_cmd2iso(df, iso,wgt)[1]

        if dist < distmin
            distmin = dist
            bestiso= iso
            bestage= age
            bestfeh= feh
            
            # save last cmd..
        end

        idiso +=1
    end

    # println("## Best distance cmd: $distmin")
    agemyr= trunc(Int, 10^bestage / 1e6)
    println("### Best age: $agemyr Myr ($bestage)")
    println("### Best FeH: $bestfeh")

    feh_gaia = median(df.mh)
    disp_feh= std(df.mh)
    println("### [Fe/H]_Gaia = $feh_gaia ($disp_feh)")

    return(bestage, bestfeh, bestiso)
end

############################################
## read the serialized MIST isochrones
## isodir: directoroy of the MIST isochrones
function read_serial_mist(isodir)
    # dirloc= pwd()
    # cd(isodir)
    arrIso=Any[]

    println("### Reading isochrone MIST models..")
    f= glob("iso_mist..*",isodir)
  
    for file in f
        df= CSV.File(file) |> DataFrame
        push!(arrIso, df)
    end
   return(arrIso)
end
###########################################
### update oc df replacing NaN with median
function update_nan_oc(df)
    println("### clean and replace missing values by median in oc solution ...")
    println("### adding new fields for magnitude..")
    nc= count(ismissing, df.ag)
    oc= filter(:ag => ag -> !any(f -> f(ag), (ismissing, isnothing, isnan)), df)
    ag0 = median(oc.ag)
    ebmr0 = median(oc.ebmr)
    mh0= median(oc.mh)
    a00= median(oc.a0)

    replace!(df.ag, NaN => ag0)
    replace!(df.ebmr, NaN => ebmr0)
    replace!(df.mh, NaN => mh0)
    replace!(df.a0, NaN => a00)

    df.gbar0 = df.gbar .- df.ag
    df.BmR0 = df.bp .- df.rp .- df.ebmr
    ## absolute magnitude
    df.G = df.gbar0 .+ 5 .- 5log10.(df.distance)


    return(df)
end
##########################################
### get solar mass from isochrone
### 
function get_star_mass(df, iso)
    global n
    println("### add star mass...")

    n= 0
    try
        # debug_red(iso)
        global n= length(iso.Gaia_G_EDR3)
        debug_red("n: $n")
    catch e
        println("### Warning : no valid Gaia_G_EDR3 field in solution...")
        return(0)
    end

    diso= zeros(2,n)
    diso[1,:] = iso.Gaia_BPmRP_EDR3
    diso[2,:] = iso.Gaia_G_EDR3

    kdtree = KDTree(diso; leafsize = 10)

    npts= length(df.G)
    pts= zeros(2,npts)
    pts[1,:] = df.BmR0
    pts[2,:] = df.G

    idxs, dists = knn(kdtree, pts, 1, true)
    idx= collect(Iterators.flatten(idxs))

    df.star_mass= iso.star_mass[idx]

    return(df)
end
################################################
## wrapper function to perform fit_isochrone fitting
## df: oc solution for the cycle
function perform_isochrone_fitting(df, isomodeldir)
    debug_red("perform_isochrone_fitting...")

    df= update_mag(df,false)
    oc= filter(:G => G -> !any(f -> f(G), (ismissing, isnothing, isnan)), df)

    noc= nrow(oc)
    debug_red("...length $noc")

    ## reading isochrone model
    if isomodeldir == "" 
        println("## Set to default isochrone models...")
        rootdir =  ENV["GAIA_ROOT"]
        isomodeldir= joinpath(rootdir,"run/data/isochrones")
    end

    debug_red(isomodeldir)
    arrIso = read_serial_mist(isomodeldir)

    if size(arrIso)[1] == 0
            println("### Isochrone model in $isomodeldir not set properly...")
            exit()
    end

    wgt= weight_cmd(oc, 0.25, 1)

    age, feh, iso= fit_isochrone(oc,arrIso, wgt)
    

    df= update_nan_oc(df)

    df= get_star_mass(df, iso)

    if typeof(df) == Int64          ## no solution...
        println("### Warning: no solution...")
        feh_gaia= -99
    else
        feh_gaia= median(df.mh)
    end

    return(df , age, feh, feh_gaia, iso)
end
## counting valid row in a df...
function counting_valid(df)
    dft= filter(:mh => mh -> !any(f -> f(mh), (ismissing, isnothing, isnan)), df)

    nvalid= length(dft.mh)

    return(nvalid)
end