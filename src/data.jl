## from module gaiaClustering
## Functions to deal with GAIA data and to normalize
##

############
## Query the GAIA data towards the coord with a radius conesearch
#####
function query_gaia(coord, radius, dump = false)
#####
    return(0)

end

function copy(s::Df)::Df
    c = Df(s.ndata, zeros(length(s.data[:,1]),s.ndata),zeros(length(s.raw[:,1]),s.ndata), zeros(length(s.err[:,1]),s.ndata))
    c.data[:,:] = s.data[:,:]
    c.raw[:,:]  = s.raw[:,:]
    c.err[:,:]  = s.err[:,:]

    return(c)
end

## dummy ...
function copy1(s::Df)::Df
    c = Df(s.ndata, zeros(length(s.data[:,1]),s.ndata),zeros(length(s.raw[:,1]),s.ndata), zeros(length(s.err[:,1]),s.ndata))
    c.data[:,:] = s.data[:,:]
    c.raw[:,:]  = s.raw[:,:]
    c.err[:,:]  = s.err[:,:]

    return(c)
end

#########
## function to create the df
######
function read_votable(voname::String)
######
    warnings = pyimport("warnings")
    warnings.filterwarnings("ignore")
    votable = pyimport("astropy.io.votable")
    vot = votable.parse(voname)
    data = vot.get_first_table()

    println("## Votable $voname read")

    ### return(data["array"]["data"])
    return(data.array.data)
end


#########
function filter_data(gaia, dist_range = [0., 2000], vra_range = [-250,250], vdec_range = [-250.,250], mag_range =[-1e9, 1e9])::Df
########
    ngaia = length(gaia)

    source_id = zeros(ngaia)
    lgal = zeros(ngaia)
    bgal = zeros(ngaia)
    ra = zeros(ngaia)
    dec= zeros(ngaia)
    distance = zeros(ngaia)
    pmra = zeros(ngaia)
    pmdec = zeros(ngaia)
    parallax = zeros(ngaia)
    vra = zeros(ngaia)
    vdec = zeros(ngaia)
    g = zeros(ngaia)
    rp = zeros(ngaia)
    bp = zeros(ngaia)
    parallax_error = zeros(ngaia)
    pmra_error     = zeros(ngaia)
    pmdec_error    = zeros(ngaia)
    radialvel      = zeros(ngaia)
    ## Galactic proper motion and velocities
    pml = zeros(ngaia)
    pmb = zeros(ngaia)
    vl = zeros(ngaia)
    vb = zeros(ngaia)

    ## Extinction A_G
    ag= zeros(ngaia)

    for i in 1:ngaia

        source_id[i]= convert(Float64, get(gaia,i-1).source_id)
        lgal[i]     = convert(Float64, get(gaia,i-1).l)
        bgal[i]     = convert(Float64, get(gaia,i-1).b)
        ra[i]       = convert(Float64, get(gaia,i-1).ra)
        dec[i]      = convert(Float64, get(gaia,i-1).dec)
        distance[i] = 1000. / convert(Float64, get(gaia,i-1).parallax)
        parallax[i] = convert(Float64, get(gaia,i-1).parallax)
        pmra[i]     = convert(Float64, get(gaia,i-1).pmra)
        pmdec[i]    = convert(Float64, get(gaia,i-1).pmdec)
        vra[i]      = 4.74e-3 * pmra[i]  * distance[i]
        vdec[i]     = 4.74e-3 * pmdec[i] * distance[i]

        ## Galactic Proper motions
        muG = PM_equatorial2galactic(pmra[i],pmdec[i]  , ra[i] , dec[i] , lgal[i])
        pml[i] = muG[1]
        pmb[i] = muG[2]
        vl[i]  = 4.74e-3 * pml[i]  * distance[i]
        vb[i]  = 4.74e-3 * pmb[i]  * distance[i]

        #fix for EDR3
        radialvel[i]    = convert(Float64, get(gaia,i-1).dr2_radial_velocity)

        ### errors.
        parallax_error[i]  = convert(Float64, get(gaia,i-1).parallax_error)
        pmra_error[i]  = convert(Float64, get(gaia,i-1).pmra_error)
        pmdec_error[i] = convert(Float64, get(gaia,i-1).pmdec_error)

        g[i]        = convert(Float64, get(gaia,i-1).phot_g_mean_mag)
        rp[i]       = convert(Float64, get(gaia,i-1).phot_rp_mean_mag)
        bp[i]       = convert(Float64, get(gaia,i-1).phot_bp_mean_mag)

        #extinction
        # fix fo EDR3, extinction not found
        ag[i]       = 99.  ##  convert(Float64, get(gaia,i-1).a_g_val)

    end

    ## Filtering ...
    i1 =  distance .> dist_range[1]
    i2 =  distance .< dist_range[2]
    i3 =  vra .> vra_range[1]
    i4 =  vra .< vra_range[2]
    i5 =  vdec .> vdec_range[1]
    i6 =  vdec .< vdec_range[2]
    i7 = g  .> mag_range[1]
    i8 = g  .< mag_range[2]
    i9  = rp  .> mag_range[1]
    i10 = rp  .< mag_range[2]
    i11  = bp  .> mag_range[1]
    i12  = bp  .< mag_range[2]

    ifinal = i1 .& i2 .& i3 .& i4 .& i5 .& i6 .& i7 .& i8 .& i9 .& i10 .& i11 .& i12

    ## G magnitude
    gbar =  g[ifinal] - 5 .* log10.(distance[ifinal]) .+ 17.

    ## Df of the filtered dat
    ndata = length(distance[ifinal])
    s = Df(ndata, zeros(8,ndata), zeros(15,ndata) , zeros(8,ndata) )

    s.data[1,:] = lgal[ifinal]
    s.data[2,:] = bgal[ifinal]
    s.data[3,:] = distance[ifinal]
    s.data[4,:] = vl[ifinal]
    s.data[5,:] = vb[ifinal]
    s.data[6,:] = gbar
    s.data[7,:] = g[ifinal] .- rp[ifinal]
    s.data[8,:] = bp[ifinal] .- g[ifinal]

    s.raw[1,:] = ra[ifinal]
    s.raw[2,:] = dec[ifinal]
    s.raw[3,:] = lgal[ifinal]
    s.raw[4,:] = bgal[ifinal]
    s.raw[5,:] = parallax[ifinal]
    s.raw[6,:] = pmra[ifinal]
    s.raw[7,:] = pmdec[ifinal]
    s.raw[8,:] = pml[ifinal]
    s.raw[9,:] = pmb[ifinal]
    s.raw[10,:] = g[ifinal]
    s.raw[11,:] = rp[ifinal]
    s.raw[12,:] = bp[ifinal]
    s.raw[13,:] = radialvel[ifinal]
    s.raw[14,:] = ag[ifinal]
    s.raw[15,:] = source_id[ifinal]

    ## Errors ..
    s.err[1,:] = parallax_error[ifinal]
    s.err[4,:] = pmra_error[ifinal]
    s.err[5,:] = pmdec_error[ifinal]

    println("## Filtering done ...")
    println("## Stars selected: $ndata")

    return(s)
end

function equatorial2galactic(α , δ)
    ## NGP coordinates
    αG = 192.85948
    δG = 27.12825
    lascend= 32.93

    b= asind(cosd(δ)*cosd(δG)*cosd(α-αG)+sind(δ)*sind(δG))

    x= cosd(δ)*cosd(δG)*sind(α-αG)
    y= sind(δ)-sind(b)*sind(δG)
    l= atand(y,x) + lascend
    println(x)
    println(y)

    if x>=0 && y<=0  l += 360.0 end
    if x<=0 && y<0   l += 360.0 end
    #if x<0  && y>=0  l -= 180.0 end

    return(l,b)
end

## angle between two points on a sphere
function angle4sphere(long1, lat1, long2, lat2)
    dLon= long2 - long1
    cosang= cosd(lat1)*cosd(lat2)*cosd(dLon)+sind(lat1)*sind(lat2)
    ang= acosd(cosang)
end

## Transform PM from equatorial to galactic system.
## See Poleski 1997 / arXiv
## PM,,corr is from Conrad (2015)
function PM_equatorial2galactic(μα , μδ , α , δ , l )
    ## NGP coordinates
    αG = 192.85948
    δG = 27.12825

    C1 = sind(δG)*cosd(δ) - cosd(δG)*sind(δ)*cosd(α - αG)
    C2 =  cosd(δG)*sind(α - αG)
    k = 1 / sqrt(C1^2 + C2^2)
    A = k * [C1 C2 ; -C2 C1]
    PMG = A * [μα ; μδ ]

    ## PM along gal. lat. corrected for differential velocity
    ## Oort constants. Not applied.
    #A , B = (14.5 , -13.) ./ 4.74
    # PMG[1] =  PMG[1] - (A*cosd(2l) + B)

    return(PMG)
end

## Compute the X,Y,Z galactic coordinates (centered on the Galactic Center)
## See Ellsworth-Bowers et al. (2013)
## Rgal was updated from Anderson et al. (2018)
## xg,yg,zg in pc

function galXYZ(l,b,distance)
    Rgal= 8.34e3
    zsun= 25
    θ= asin(zsun/Rgal)

    xg= Rgal*cos(θ)-distance*(cosd(l)cosd(b)cos(θ)+sind(b)sin(θ))
    yg= -distance*sind(l)cosd(b)
    zg= Rgal*sin(θ)-distance*(cosd(l)cosd(b)sin(θ)-sind(b)cos(θ))

    return(xg,yg,zg)
end


### Correction of the radial velocity
### Conrad (3025)

function RVEL_corr(rvel , distance , l)
  ## Oort's constant
  A = 14.5
  rv = rvel - A * 1e-3 * distance * sind(2l)
  return(rv)
end

### compute galactic U V W from idlastro gal_uvw.pro
##
## ra, dec: degrees
## distance; pc
## pmra,pmdec: milliarcsec/yr
## vrad: km/s

## UVW: km/s
#      U - Velocity (km/s) positive toward the Galactic *anti*center
#      V - Velocity (km/s) positive in the direction of Galactic rotation
#      W - Velocity (km/s) positive toward the North Galactic Pole
#
function galUVW(ra, dec, distance, pmra, pmdec, vrad ; LSR_vel=[-8.5 ; 13.38 ; 6.49])
    k = 4.74047     #Equivalent of 1 A.U/yr in km/s

    T=  [0.0548756   0.873437  0.483835 ;
          0.494109   -0.44483   0.746982 ;
         -0.867666   -0.198076  0.455984]

    A1= [ cosd(ra) sind(ra) 0 ; sind(ra) -cosd(ra) 0 ; 0 0 -1]
    A2= [ cosd(dec) 0 -sind(dec) ; 0 -1 0 ; -sind(dec) 0 -cosd(dec)]
    vec1 = vrad
    vec2 = k*pmra*1e-3*distance
    vec3 = k*pmdec*1e-3*distance
    v= [vec1 ; vec2 ; vec3]

    uvw= T*A1*A2*v + LSR_vel

    return(uvw)
end


######
function add_cartesian(s::Df, centering = true)::Df
######
    dfresult = copy(s)
    off = zeros(2)

    if centering
        off[1] = mean(s.data[1,:])
        off[2] = mean(s.data[2,:])
    end

    lgal = DEG2RAD .* (s.data[1,:] .- off[1])
    bgal = DEG2RAD .* (s.data[2,:] .- off[2])

    dfresult.data[1,:] = s.data[3,:] .* cos.(bgal) .* cos.(lgal)
    dfresult.data[2,:] = s.data[3,:] .* cos.(bgal) .* sin.(lgal)
    dfresult.data[3,:] = s.data[3,:] .* sin.(bgal)

    println("## Cartesian transformation done ...")

    return(dfresult)
end


######
function  normalization_PerBlock(s::Df, block , weightblock, norm , density = false , verbose = true)
######
    dfresult = copy(s)
    ndf = size(s.data)
    scale8d = zeros(ndf[1])
    vector8d = 0.

    ind = 1
    for aw in zip(block,weightblock)
        weight = aw[2]
        for ak in aw[1]
            normK = normalizationVector(norm, density, dfresult.data[ak,:])
            # normK[2] = normK[2] * totalWeight
            dfresult.data[ak,:]    =   weight .* (s.data[ak,:] .- normK[1] ) ./ normK[2]
            scale8d[ind] = weight / normK[2]
            vector8d += scale8d[ind] ^ 2
            ind += 1
        end
    end

    vector8d = sqrt(vector8d)
    scale8d[:] = scale8d[:] ./ vector8d
    dfresult.data[:,:] = dfresult.data[:,:] ./ vector8d

    if verbose
        println("## Normalization $norm done...")
        println("### [1pc,1pc,1pc,1km/s,1km/s,1mag,1mag,1mag] equivalent to $scale8d")
        println("##")
    end

    return(dfresult , scale8d)
end

######
function normalizationVector(norm, density, arr)
######
    vecNorm = [0.0, 1.0]

    if norm == "identity"
        vecNorm = [0.0, 1.0]

        elseif norm == "normal"
            stdArr = std(arr)
            meanArr = mean(arr)
            vecNorm  = [meanArr , stdArr]

        elseif norm == "minmax"
            minarr  = minimum(vcat(arr...))
            maxarr  = maximum(vcat(arr...))
            vecNorm = [minarr, maxarr-minarr]

    end

    if density
        vecNorm[2] = vecNorm[2] * length(arr)
    end

    return(vecNorm)
end

######
function subsetDf(df::Df, indx)::Df
######

  ndat = length(indx)
  subset = Df(ndat, df.data[:,indx], df.raw[:,indx], df.err[:,indx])

    return(subset)
end


####
# Create the DataFrame to save the cluster...
##
function export_df(votname, ocdir, df , dfcart, labels , labelmax)
    ra= df.raw[1, labels[labelmax]]
    dec= df.raw[2,labels[labelmax]]
    l= df.data[1, labels[labelmax]]
    b= df.data[2,labels[labelmax]]
    d= df.data[3,labels[labelmax]]
    pmra= df.raw[6, labels[labelmax]]
    pmdec= df.raw[7, labels[labelmax]]
    X= dfcart.data[1, labels[labelmax]]
    Y= dfcart.data[2, labels[labelmax]]
    Z= dfcart.data[3, labels[labelmax]]
    vl= df.data[4,labels[labelmax]]
    vb= df.data[5,labels[labelmax]]
    vrad= df.raw[13,labels[labelmax]]
    gbar= df.raw[10,labels[labelmax]]
    rp= df.raw[11,labels[labelmax]]
    bp= df.raw[12,labels[labelmax]]
    ag= df.raw[14,labels[labelmax]]

    source_id= df.raw[15,labels[labelmax]]

    oc= DataFrame(sourceid=source_id,ra=ra,dec=dec,l=l,b=b, distance=d,pmra=pmra, pmdec=pmdec, X=X,Y=Y,Z=Z,vl=vl,vb=vb,vrad=vrad,gbar=gbar,rp=rp,bp=bp, ag=ag)

    name= split(votname,".")
    infix= ""
    for iname in name
        if iname != "vot"
            infix *= iname*"."
        end
    end
    infix *= "oc.csv"
    filename= @sprintf("%s/%s",ocdir, infix)
    CSV.write(filename,oc,delim=';')
    @printf("### %s created \n",filename)
end
