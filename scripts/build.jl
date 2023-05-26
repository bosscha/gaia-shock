## standalone script to (re)build sets of candidate cluster (cycle) from different
## fields
##

using DataFrames , CSV, TOML , ArgParse
using Dates, Printf, PyCall

rootdir =  ENV["GAIA_ROOT"]

push!(LOAD_PATH,"$rootdir/run/src")
using GaiaClustering

#########################  reprocess function
function reprocess(meta)

    tstart= now()
    println(blue("## Reprocessing ..."))
    println(blue("## Starting at $tstart"))

    mrepro= meta["reprocess"]
    mgene= meta["general"]
    mextra= read_params(mrepro["extrafile"], false)

    
    checkdir(joinpath(mgene["wdir"],mextra.ocdir), joinpath(mgene["wdir"],mextra.plotdir))
    
    debug_red(joinpath(mgene["wdir"],mextra.plotdir))

    cd(mgene["wdir"])

    progressfile= "_done.csv"            #progress file 

    # if !haskey(mextra, :rootdir) mextra.rootdir= "./" end
    # if !haskey(mextra, :wdir) mextra.wdir= "./" end
    #if !haskey(mextra, :plotdir) mextra.plotdir= "./plotSelect" end
    # if !haskey(mextra, :ocdir) mextra.ocdir= "./oc" end

    println(mrepro)
    println(mgene)

    if isfile(progressfile)              
        dfp=  CSV.File(progressfile, delim=",") |> DataFrame 
    else
        dfp= DataFrame(votname= String[])
    end

    dfblck= get_blacklist(mgene)

    Kdeg= 57.69

    if mextra.optim == "no"
        optim= false
        println("## Using preprocessed file : $(mrepro["optsol"])")
        ## delim Guess
        fline= readline(mrepro["optsol"]) 
        if occursin(";",fline) delimiter= ";" else delimiter= "," end
        dfoptsol= CSV.File(mrepro["optsol"], delim=delimiter) |> DataFrame
        
        gaia= pyimport("astroquery.gaia")

        radone= [] ; decdone= []
        for row in eachrow(dfoptsol)
            if mgene["getvot"] == "yes"
                rect= false
                ra= row.ra ; dec= row.dec
                push!(radone,ra) ; push!(decdone, dec)
                tol= mrepro["tol"] ; radius= mrepro["radius"]

                estangle= Kdeg * 25 / row.distance             # estimated field size in degree (if angle small for tangent)
                debug_red("est. angle $estangle deg")

                name= @sprintf("RA%.3fDec%.3f",ra,dec)
                votname= @sprintf("%s-%2.1fdeg.vot",name, radius)

                if votname in dfp.votname || votname in dfblck.votname
                    println("## $votname skipped..")
                else

                    mextra.votname= get_gaia_data_many(gaia, radius, tol, ra, dec, name , rect)

                    if haskey(row, :w3dm)    mextra.w3d= row.w3dm else mextra.w3d= row.w3d end
                    if haskey(row, :wvelm)   mextra.wvel= row.wvelm else mextra.wvel= row.wvel end     
                    if haskey(row, :whrdm)   mextra.whrd= row.whrdm else mextra.whrd= row.whrd end
                    if haskey(row, :epsm)    mextra.eps= row.epsm else mextra.eps= row.eps end
                    if haskey(row, :mclm)    mextra.mcl= Int(floor(row.mclm)) else mextra.mcl= Int(floor(row.mcl)) end
                    if haskey(row, :mneim)    mextra.mnei= Int(floor(row.mneim)) else mextra.mnei= Int(floor(row.mnei)) end
               
                    debug_red(mextra)
                    extra(mextra, optim)
                    push!(dfp, [mextra.votname])
                    CSV.write(progressfile, dfp, delim=";")

                    if mgene["rmvot"] == "yes"
                        rm(mextra.votname)
                        println("## votable $(mextra.votname) removed")
                    end
                end
                plot_sky(radone, decdone, radius=50, figname= "reprocess-allsky.png")
            end
        end
    end
    rm(progressfile)
end

#########################  randomfields function
function randomfields(meta)
    tstart= now()
    println(blue("## Processing random fields ..."))
    println(blue("## Starting at $tstart"))

    mrandom= meta["random"]
    mgene= meta["general"]
    mextra= read_params(mrandom["extrafile"], false)

    cd(mgene["wdir"])
    progressfile= "_done.csv"            #progress file 
    mextra.rootdir= "./"
    mextra.wdir= "./"
    mextra.plotdir= "./plotSelect"
    mextra.ocdir= "./oc"

    println(mrandom)
    println(mgene)

    if isfile(progressfile)              
        dfp=  CSV.File(progressfile, delim=",") |> DataFrame 
    else
        dfp= DataFrame(votname= String[])
    end

    ndone= size(dfp)[1]
    nfields= mrandom["fields"]
    radius= mrandom["radius"]
    tol= mrandom["tol"]
    mode= mrandom["mode"]
    
    if mode == "galactic"
        bscale= mrandom["bscale"]
    else
        bscale= -1
    end
    rect= false   # conesearch

    if ndone < nfields
        notfinished= true
    else
        notfinished= false
    end

    gaia= pyimport("astroquery.gaia")

    if mextra.optim == "yes"
        optim= true
    elseif mextra.optim == "no"
        optim= false
    end

    while notfinished
        ra, dec = get_random_field(mode, bscale)
        name= @sprintf("RA%.3fDec%.3f",ra,dec)
        debug_red(name)

        mextra.votname= get_gaia_data_many(gaia, radius, tol, ra, dec, name , rect)
        extra(mextra, optim)

        push!(dfp, [mextra.votname])
        CSV.write(progressfile, dfp, delim=";")
        if mgene["rmvot"] == "yes"
            rm(mextra.votname)
            println("## votable $(mextra.votname) removed")
        end

        ndone += 1
        println("## $ndone random fields processed...")
        if ndone > nfields
            notfinished= false
        end    
    end
    rm(progressfile)
end

#########################  gridding function
function gridding(meta)

    tstart= now()
    println(blue("## Processing gridding fields ..."))
    println(blue("## Starting at $tstart"))

    mgrid= meta["gridding"]
    mgene= meta["general"]
    mextra= read_params(mgrid["extrafile"], false)

    cd(mgene["wdir"])
    progressfile= "_done.csv"            #progress file 
    mextra.rootdir= "./"
    mextra.wdir= "./"
    mextra.plotdir= "./plotSelect"
    mextra.ocdir= "./oc"

    println(mgrid)
    println(mgene)

    if isfile(progressfile)              
        dfp=  CSV.File(progressfile, delim=",") |> DataFrame 
    else
        dfp= DataFrame(votname= String[])
    end

    gaia= pyimport("astroquery.gaia")

    radius= mgrid["radius"]
    tol= mgrid["tol"]
    ref= mgrid["ref"]
    rect= false   # conesearch

    if ref == "galactic"
        xiter= mgrid["lgal_grid"]
        yiter= mgrid["bgal_grid"]
    elseif ref == "equatorial"
        xiter= mgrid["ra_grid"]
        yiter= mgrid["dec_grid"]  
    end
    
    if mextra.optim == "yes"
        optim= true
    elseif mextra.optim == "no"
        optim= false
    end
    
    debug_red(xiter)
    radone= [] ; decdone= []

    for xx in xiter[1]:xiter[3]:xiter[2]
        for yy in yiter[1]:yiter[3]:yiter[2]
            if ref == "galactic"
                println("## Starting (l,b): $xx $yy")
                ra, dec= galactic2equatorial(xx,yy)
            elseif ref == "equatorial"
                println("## Starting (RA,Dec): $xx $yy")
                ra= xx ; dec= yy
            end
            push!(radone,ra)
            push!(decdone,dec)           

            name= @sprintf("RA%.3fDec%.3f",ra,dec)
            votname= @sprintf("%s-%2.1fdeg.vot",name, radius)

            if votname in dfp.votname 
                println("## $name done...")
            else
                mextra.votname= get_gaia_data_many(gaia, radius, tol, ra, dec, name , rect)
                extra(mextra, optim)
    
                push!(dfp, [mextra.votname])
                CSV.write(progressfile, dfp, delim=";")
                if mgene["rmvot"] == "yes"
                    rm(mextra.votname)
                    println("## votable $(mextra.votname) removed")
                end
            end           
            plot_sky(radone,decdone , radius=20, figname="gridding-allsky.png")
        end
    end
end
#########################  merge function
function merge(meta)

    tstart= now()
    println(blue("## Merging catalog ..."))
    println(blue("## Starting at $tstart"))

    mmerge= meta["merge"]
    mgene= meta["general"]

    cd(mgene["wdir"])

    catalog= mmerge["catalog"]
    mergefile= name= @sprintf("%s.merge", catalog)
    mode= mmerge["mode"]
    
    debug_red(mergefile)

    if mode == "duplicate"
        println("### Merge, removing duplicated clusters...")
        debug_red(mmerge)

        toldeg= mmerge["toldeg"]
        toldist= mmerge["toldist"]
        tolndiff= mmerge["tolndiff"]
        metric= mmerge["metric"] 
        println("### Merge, metric $metric")
        
        dfcat=  CSV.File(catalog, delim=";") |> DataFrame

        dfmerge= rm_duplicated(dfcat, toldeg, toldist, tolndiff, metric)
        CSV.write(mergefile, dfmerge, delim=";")
        println("## Catalog $mergefile created.")
    end 

    if mode == "simbad"
        println("### Merge, searching for objects with Simbad...")
        coord= pyimport("astropy.coordinates")
        Simbad= pyimport("astroquery.simbad")

        tolquery= mmerge["tolquery"]       ## arcmin radius for the object search

        dfcat=  CSV.File(catalog, delim=";") |> DataFrame

        namecla= []
        for row in eachrow(dfcat)
            ra= row.ra ; dec= row.dec
            c= coord.SkyCoord(ra,dec, unit="deg")
            res= Simbad.Simbad.query_region(c, radius="$tolquery arcmin")
            namecl= "-"

            if res !== nothing
                for lobj in res
                    t= split(lobj[1])
                    if t[1]=="Cl*"
                        namecl= t[2]*"_"*t[3]
                        debug_red(namecl)
                    end
                end
            end
            push!(namecla,namecl)
        end
        dfcat[!,:name] = namecla

        CSV.write(mergefile, dfcat, delim=";")
        println("## Catalog $mergefile created.")
    end

end
#########################
function get_blacklist(m)
    if haskey(m, "blacklist")
        blackfile= m["blacklist"]
        if isfile(blackfile)              
            dfblck=  CSV.File(blackfile, delim=",") |> DataFrame 
            println("### Blacklist $blackfile read")
        else
            dfblck= DataFrame(votname= [""])
        end
    else
        dfblck= DataFrame(votname= [""])
    end
    return(dfblck)
end
#################################### MAIN ########################
let
    println(ARGS)
    println("############################")
    println("### Building Gaia results...")

    metabuild = TOML.parsefile(ARGS[1])
    
    key = collect(keys(metabuild))
    println(key)
    for k in key
        # key1 = collect(keys(metabuild[k]))

        if k == "reprocess" 
            reprocess(metabuild)
        elseif k == "random" 
            randomfields(metabuild)
        elseif k== "gridding"
            gridding(metabuild)
        elseif k== "merge"
            merge(metabuild)
        end       
    end
end

