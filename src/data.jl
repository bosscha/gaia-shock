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
    warnings[:filterwarnings]("ignore")
    votable = pyimport("astropy.io.votable")
    vot = votable[:parse](voname)
    data = vot[:get_first_table]() 
    
    println("## Votable $voname read")
    
    return(data["array"]["data"])
end


#########
function filter_data(gaia, dist_range = [0., 2000], vra_range = [-250,250], vdec_range = [-250.,250], mag_range =[-1e9, 1e9])::Df
########
    ngaia = length(gaia)
    
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
    
    for i in 1:ngaia
        lgal[i]     = convert(Float64,gaia[i]["l"])
        bgal[i]     = convert(Float64,gaia[i]["b"])
        ra[i]       = convert(Float64,gaia[i]["ra"])
        dec[i]      = convert(Float64,gaia[i]["dec"])
        distance[i] = 1000. / convert(Float64,gaia[i]["parallax"])
        parallax[i] = convert(Float64,gaia[i]["parallax"])
        pmra[i]     = convert(Float64,gaia[i]["pmra"])
        pmdec[i]    = convert(Float64,gaia[i]["pmdec"])
        vra[i]      = 4.74e-3 * pmra[i]  * distance[i]
        vdec[i]     = 4.74e-3 * pmdec[i] * distance[i]
        ### errors.
        parallax_error[i]  = convert(Float64,gaia[i]["parallax_error"])
        pmra_error[i]  = convert(Float64,gaia[i]["pmra_error"])
        pmdec_error[i] = convert(Float64,gaia[i]["pmdec_error"])
        
        g[i]        = convert(Float64,gaia[i]["phot_g_mean_mag"])
        rp[i]       = convert(Float64,gaia[i]["phot_rp_mean_mag"])
        bp[i]       = convert(Float64,gaia[i]["phot_bp_mean_mag"])
        
        
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
     # gbar =  g[ifinal] - 5. * log10.(distance[ifinal]) + 17.
    gbar =  g[ifinal] - 5 .* log10.(distance[ifinal]) .+ 17.
    
    ## Df of the filtered dat
    ndata = length(distance[ifinal])
    s = Df(ndata, zeros(8,ndata), zeros(10,ndata) , zeros(8,ndata) )
    
    s.data[1,:] = lgal[ifinal]
    s.data[2,:] = bgal[ifinal]
    s.data[3,:] = distance[ifinal]
    s.data[4,:] = vra[ifinal]
    s.data[5,:] = vdec[ifinal]
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
    s.raw[8,:] = gbar
    s.raw[9,:] = rp[ifinal]
    s.raw[10,:] = bp[ifinal] 
 
    
    ## Errors ..
    s.err[1,:] = parallax_error[ifinal]
    s.err[4,:] = pmra_error[ifinal]
    s.err[5,:] = pmdec_error[ifinal]
    
    println("## Filtering done ...")
    println("## Stars selected: $ndata")
    
    return(s)
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


