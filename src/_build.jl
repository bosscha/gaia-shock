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