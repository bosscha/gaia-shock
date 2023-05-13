## functions to be used  by the build reprocessing scripts

function extra(m::GaiaClustering.meta, optim)
    tstart= now()
    rng = MersenneTwister()
    uuid=uuid4(rng)
    m.uuid= uuid

    println("###########################")
    println("## Starting with $(m.votname)")
    println("## Starting at $tstart")
    println("## Id $uuid")

    df , dfcart , dfcartnorm = get_data(m)

    cycle, flag= cycle_extraction_optim(df, dfcart, m, optim)

    tend= now()
    println("## Ending at $tend")
    println("## number of cycle: $cycle , flag:$flag ")
    println("##########################")
    println("##")

    duration= Dates.value(tend-tstart) / (1000*3600)
    durationstr= @sprintf("%3.3f", duration)
    @printf("## %s \n",specialstr("Duration: $durationstr hours","YELLOW"))
    @printf("## %s \n",specialstr("Votable done: $(m.votname)","YELLOW"))
    println("##\n##")

    debug_red("Cleaning ...")
    df=0 ; dfcart= 0 ; dfcartnorm= 0
    GC.gc()
end
########
function get_gaia_data(radius, tol, ra, dec, name, rect, table= "gaiadr3.gaia_source")

    gaia= pyimport("astroquery.gaia")
    # u= pyimport("astropy.units")
    # coord= pyimport("astropy.coordinates")

    if rect
        adql= @sprintf("SELECT * FROM %s WHERE %s.ra BETWEEN %f AND %f AND
%s.dec BETWEEN %f AND %f AND abs(pmra_error/pmra)<  %3.3f AND abs(pmdec_error/pmdec)< %3.3f
    AND abs(parallax_error/parallax)< %3.3f;", table, table, ra-radius, ra+radius, table,dec-radius,
    dec+radius,tol, tol, tol )
    else
        adql= @sprintf("SELECT * FROM %s WHERE CONTAINS(POINT('ICRS',%s.ra, %s.dec),
    CIRCLE('ICRS',%f, %f, %f)) = 1  AND abs(pmra_error/pmra)<  %3.3f AND abs(pmdec_error/pmdec)< %3.3f
    AND abs(parallax_error/parallax)< %3.3f;", table, table, table, ra, dec, radius, tol, tol, tol )
    end

    println("## Downloading data...")
    println(adql)
    job= gaia.Gaia.launch_job_async(adql, dump_to_file=true)
    outfile= job.outputFile
    filedst= @sprintf("%s-%2.1fdeg.vot",name, radius)
    mv(outfile, filedst , force=true)
    println("## Votable saved in: $filedst")
    
    return(filedst)
end

########
## call many times to the gaia archive. gaia object passed
function get_gaia_data_many(gaia, radius, tol, ra, dec, name, rect, table= "gaiadr3.gaia_source")

    if rect
        adql= @sprintf("SELECT * FROM %s WHERE %s.ra BETWEEN %f AND %f AND
%s.dec BETWEEN %f AND %f AND abs(pmra_error/pmra)<  %3.3f AND abs(pmdec_error/pmdec)< %3.3f
    AND abs(parallax_error/parallax)< %3.3f;", table, table, ra-radius, ra+radius, table,dec-radius,
    dec+radius,tol, tol, tol )
    else
        adql= @sprintf("SELECT * FROM %s WHERE CONTAINS(POINT('ICRS',%s.ra, %s.dec),
    CIRCLE('ICRS',%f, %f, %f)) = 1  AND abs(pmra_error/pmra)<  %3.3f AND abs(pmdec_error/pmdec)< %3.3f
    AND abs(parallax_error/parallax)< %3.3f;", table, table, table, ra, dec, radius, tol, tol, tol )
    end

    println("## Downloading data...")
    println(adql)
    job= gaia.Gaia.launch_job_async(adql, dump_to_file=true)
    outfile= job.outputFile
    filedst= @sprintf("%s-%2.1fdeg.vot",name, radius)
    mv(outfile, filedst , force=true)
    println("## Votable saved in: $filedst")
    
    return(filedst)
end

###########
### create random field positions
### mode:uniform|galactic   if galactic needs to set b scale (hb)
###
function get_random_field(mode="uniform", hb=0)
    if mode=="uniform"
        ra= rand(Float64)*360.
        dec= rand(Float64)*180 - 90.
    elseif mode=="galactic"
        coord= pyimport("astropy.coordinates")
        l= rand(Float64)*360 -180
        b= hb*randexp()
        if rand() < 0.5 b= -b end
        println("### Random field: l:$l b:$b ")
        radec= coord.SkyCoord(l,b, unit="deg", frame="galactic").icrs
        ra= radec.ra[1] ; dec= radec.dec[1]
    else
        println("## Warning! Random field mode not found")
        return(0,0)
    end
    return(ra,dec)
end
##############
function galactic2equatorial(l,b)
    coord= pyimport("astropy.coordinates")
    radec= coord.SkyCoord(l,b, unit="deg", frame="galactic").icrs
    ra= radec.ra[1] ; dec= radec.dec[1]
    return(ra,dec)
end
################
## remove duplicated oc from a df
## 
## df: catalog DataFrame
## toldeg: tolerance in degree for the position RA,dec
## toldist: tolerance in parsec for the distance
## metric: metric (Qn|Edgm) to choose the oc
## rmfile: remove permanently the oc and plot files
## tolndiff: minimum difference relative between N for two candidates
##
function rm_duplicated(df, toldeg, toldist, tolndiff, metric= "Qn", rmfile= false)
    dfmerge= df

    s= size(dfmerge)
    ndrop= []

    for i in 1:s[1]
        for j in i+1:s[1]
            ra1= dfmerge[i,"ra"] ; dec1= dfmerge[i,"dec"]
            ra2= dfmerge[j,"ra"] ; dec2= dfmerge[j,"dec"]
            dist1= dfmerge[i,"distance"]
            dist2= dfmerge[j,"distance"] 
            n1= dfmerge[i,"nstars"]
            n2= dfmerge[j,"nstars"]
            edg1= dfmerge[i,"edgratm"]
            edg2= dfmerge[j,"edgratm"]
            name1= dfmerge[i,"votname"]
            name2= dfmerge[j,"votname"]
            
            if abs(ra1-ra2) < toldeg && abs(dec1-dec2) < toldeg && abs(dist1-dist2) < toldist && name1 != name2
                if min(n1,n2) / max(n1,n2) < tolndiff
                    println("## Warning merge, two candidates have large N difference...")
                end

                if metric == "Qn"
                    if n1 > n2 && j ∉ ndrop
                        push!(ndrop,j)
                    elseif i ∉ ndrop
                        push!(ndrop,i) 
                    end
                elseif metric == "Edgm"
                    if edg2 > edg1 && j ∉ ndrop
                        push!(ndrop,j)
                    elseif i ∉ ndrop
                        push!(ndrop,i) 
                    end
                end
            end
        end
    end
    sort!(ndrop)
   
    if rmfile
        for i in ndrop
            name1= dfmerge[i,"votname"]
            cycle1= dfmerge[i,"cycle"]
            
            wildoc= @sprintf("./oc/%s.%d.*.csv",name1,cycle1)
            wildplot= @sprintf("./plotSelect/%s.%d.*.png",name1,cycle1)

            rm.(glob(wildoc))            
            rm.(glob(wildplot))
        end
        println("### Files with duplicated solutions removed.")
    end
    delete!(dfmerge,ndrop)

    sdrop= size(ndrop)[1]
    smerge= size(dfmerge)[1]
    println("### $sdrop duplicated solutions dropped.")
    println("### $smerge solutions left. ")

    return(dfmerge)
end
#######################################
##### distance of two CMDs
#### df1, df2: two solutions
function distance_cmd(df1, df2)
    debug_red(df1["gbar"])

    n= length(iso.Gaia_G_EDR3)

    dfref= zeros(2,n)
    dfref[1,:] = iso.Gaia_BPmRP_EDR3
    dref[2,:] = iso.Gaia_G_EDR3

    kdtree = KDTree(dref; leafsize = 10)

    npts= length(df.G)
    pts= zeros(2,npts)
    pts[1,:] = df.BmR0
    pts[2,:] = df.G

    idxs, dists = nn(kdtree, pts)
    total= sum(dists)
end