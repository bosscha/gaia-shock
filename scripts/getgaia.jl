## wrapper and line command script to get GAIA votable in a conesearch
using Printf, ArgParse
using PyCall

gaia= pyimport("astroquery.gaia")
u= pyimport("astropy.units")
coord= pyimport("astropy.coordinates")

##
## parse getgaia options/flags
##
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-f"
            help = "Field name"
            arg_type = String
            default= "gaia"
        "--resolv", "-r"
            help = "To resolv the field name on the sky"
            action = :store_true
        "--rad"
            help = "Radius on sky in degrees."
            arg_type = Float64
            default= 2.0
        "--tol"
            help = "Tolerance on Gaia data."
            arg_type = Float64
            default= 0.2
        "-l"
            help = "Galactic longitud"
            arg_type = Float64
        "-b"
            help = "Galactic latitude"
            arg_type = Float64
        "--ra"
            help = "Right Ascension in degrees"
            arg_type = Float64
        "--dec"
            help = "Declination in degrees"
            arg_type = Float64

    end

    return parse_args(s)
end
## query the gaia data and save then
function get_gaia_data(radius, tol, ra, dec, name, table= "gaiaedr3.gaia_source")

    adql= @sprintf("SELECT * FROM %s WHERE CONTAINS(POINT('ICRS',%s.ra, %s.dec),
    CIRCLE('ICRS',%f, %f, %f)) = 1  AND abs(pmra_error/pmra)<  %3.3f AND abs(pmdec_error/pmdec)< %3.3f
    AND abs(parallax_error/parallax)< %3.3f;", table, table, table, ra, dec, radius, tol, tol, tol )


    println("## Downloading data...")
    println(adql)
    job= gaia.Gaia.launch_job_async(adql, dump_to_file=true)
    outfile= job.outputFile
    filedst= @sprintf("%s-%2.1fdeg.vot",name, radius)
    mv(outfile, filedst , force=true)
    println("## Votable saved in: $filedst")

end
function coord_galactic(longi, lati)
    c = coord.SkyCoord(l=longi*u.degree, b=lati*u.degree,frame="galactic")
    ct= c.transform_to("fk5")
    return(ct.ra[1], ct.dec[1])
end
#################################### MAIN###############################
let
    args = parse_commandline()

    fname= args["f"]
    l= args["l"]
    b= args["b"]
    ra= args["ra"]
    dec= args["dec"]
    rad= args["rad"]
    tol= args[ "tol"]
    resolv= args[ "resolv"]

    iscoord= false
    if l != nothing && b != nothing
        longi, lati= coord_galactic(l, b)
        iscoord= true
    elseif ra != nothing && dec != nothing
        longi= ra ; lati= dec
        iscoord= true
    elseif resolv
        println("## Resolving name...")
        try
            c= coord.get_icrs_coordinates(fname)
            longi= c.ra[1] ; lati= c.dec[1]
            iscoord= true
        catch err
            @printf("## failed to resolve %s \n", fname)
        end
    end

    if iscoord
        get_gaia_data(rad, tol, longi, lati, fname)
    else
        println("No coordinates defined...")
    end
end
