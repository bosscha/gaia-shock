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
## call many time to the gaia archive. gaia object passed
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