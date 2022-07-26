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

    cd(mgene["wdir"])
    progressfile= "_done.csv"            #progress file 
    mextra.rootdir= "./"
    mextra.wdir= "./"
    mextra.plotdir= "./plotSelect"
    mextra.ocdir= "./oc"

    println(mrepro)
    println(mgene)

    if isfile(progressfile)              
        dfp=  CSV.File(progressfile, delim=",") |> DataFrame 
    else
        dfp= DataFrame(votname= String[])
    end
        
    Kdeg= 57.69

    if mextra.optim == "no"
        optim= false
        println(mrepro["optsol"])
        dfoptsol= CSV.File(mrepro["optsol"], delim=",") |> DataFrame
        # println(dfoptsol)
        
        gaia= pyimport("astroquery.gaia")

        for row in eachrow(dfoptsol)
            if mgene["getvot"] == "yes"
                rect= false
                ra= row.ra ; dec= row.dec
                tol= mrepro["tol"] ; radius= mrepro["radius"]

                estangle= Kdeg * 25 / row.distance             # estimated field size in degree (if angle small for tangent)
                debug_red("est. angle $estangle deg")

                name= @sprintf("RA%.3fDec%.3f",ra,dec)
                votname= @sprintf("%s-%2.1fdeg.vot",name, radius)

                if votname in dfp.votname
                    println("## $votname skipped..")
                else

                    mextra.votname= get_gaia_data_many(gaia, radius, tol, ra, dec, name , rect)

                    mextra.w3d= row.w3dm
                    mextra.wvel= row.wvelm
                    mextra.whrd= row.whrdm
                    mextra.eps= row.epsm
                    mextra.mcl= row.mclm
                    mextra.mnei= row.mneim
                
                    extra(mextra, optim)
                    push!(dfp, [mextra.votname])
                    CSV.write(progressfile, dfp, delim=";")

                    if mgene["rmvot"] == "yes"
                        rm(mextra.votname)
                        println("## votable $(mextra.votname) removed")
                    end
                end
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

    for xx in xiter[1]:xiter[3]:xiter[2]
        for yy in yiter[1]:yiter[3]:yiter[2]
            debug_red(xx)
            if ref == "galactic"
                ra, dec= galactic2equatorial(xx,yy)
            elseif ref == "equatorial"
                ra= xx ; dec= yy
            end

            name= @sprintf("RA%.3fDec%.3f",ra,dec)
            votname= @sprintf("%s-%2.1fdeg.vot",name, radius)
            debug_red(votname)

            if name in dfp.votname 
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
        end
    end
end
#################################### MAIN ########################
let
    println(ARGS)
    println("####################")
    println("### testing build...")

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
        end       
    end

    
    # extra(0,false)
end

