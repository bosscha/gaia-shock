## Mass segregation computation in stellar clusters.
## References : Cartwright et al. 2004, 2009 ; Allison et al. 2009 ; Olczak et al. 2011 ; Yu et al. 2011
##

## Create a graph with the selected stars
## ind: subset to select the x,y
## dtyp= "2d" or "3d"
##

##!! WARNING !! 2d is for the projection YZ since X is the distance with GAIA
##
function mst_graph(ind,x,y,z=[] , dtyp="2d")
    nxy= length(ind)
    if dtyp == "2d"
        A= Array{Float64}(undef,2,nxy)
    elseif dtyp == "3d"
       A= Array{Float64}(undef,3,nxy)
    end

    for i in 1:nxy
        if dtyp == "2d"
            A[:,i]= [y[ind[i]] z[ind[i]]]
        elseif dtyp == "3d"
            A[:,i]= [x[ind[i]] y[ind[i]] z[ind[i]]]
        end
    end

    d= Euclidean()
    p= pairwise(d, A, A, dims=2)


    src= Array{Int}(undef,0)
    dst= Array{Int}(undef,0)
    wgt= Array{Float64}(undef,0)

    for i in 1:nxy
        for j in 1:nxy
            if i != j
                push!(src,i)
                push!(dst,j)
                push!(wgt,p[i,j])
            end
        end
    end

    G= SimpleWeightedGraph(src, dst, wgt)
    kr= kruskal_mst(G, G.weights)

    return(kr)
end

## Computation of the λMST using either arithmetic or geometric average (Olczak et al. 2011)
##
## dtype= "ari" or "geo"
function lambda_mst(edges, dtyp="ari")
    nedg= length(edges)
    ltot= 0

    for w in edges
        if dtyp== "ari"
            ltot+= w.weight
        elseif dtyp== "geo"
            ltot+= log(w.weight)
        end
    end

    if dtyp== "ari"
        λmst= ltot/nedg
    elseif dtyp== "geo"
        λmst= exp(ltot/nedg)
    end
    return(λmst)
end

## estimate the sample size for the stats given n and N for the massive and total sample
## p= 99%
## Olczak et al. 2011)
function sample_size(n, N, p=0.99)
    k= log(1-p)/log(1-n/N)
    return(convert(Int,ceil(k)))
end

## Compute the mass segregation metric comparing two population set (2nd is random)
##
## nrandom: number of random set (default=1)
##
function kappa_ms(indMass, x, y , z, nrandom=1, dtyp="2d", daver= "geo")
    npop= length(indMass)
    ntot= length(x)

    ## Massive stars
    krMass= mst_graph(indMass, x,y,z , dtyp)
    λMass= lambda_mst(krMass,daver)
    ## Random set
    λArr= []
    κArr= []
    for i in 1:nrandom
        perm= randperm(ntot)
        indRef= perm[1:npop]
        krRef= mst_graph(indRef, x,y,z, dtyp)
        λRef= lambda_mst(krRef,daver)
        push!(λArr,λRef)
        push!(κArr, log(λRef/λMass))
    end
    λMean= mean(λArr)
    σκ= std(κArr)

    ## ratio (log) of the reference over massive ones
    κms= log(λMean/λMass)

    return(κms , σκ)
end


## Select the N percentile most massive stars (brighter)
##
## reverse:true -> light stars
function select_massivestars(mag, percentile=15 , reverse=false)
    nstar= length(mag)
    isort= sortperm(mag, rev=reverse)
    nselect= convert(Int, floor(nstar*percentile/100))

    return(isort[1:nselect])
end


## return the MS kappa and kappa_err
## IN
## mag: magnitude array
## X,Y,Z: positions array
## perc: percentile of massive stars
## dtyp: "3d" or "2d"
## daver: "geo" or "ari"
## reverse: false:massive, true:light - stars
##
function get_kappaMS(mag, X, Y, Z, perc, dtyp="3d", daver="geo", reverse=false)
    imass= select_massivestars(mag, perc, reverse)
    ns= sample_size(length(imass), length(X))
    κms, κms_err= kappa_ms(imass, X, Y, Z, ns, dtyp, daver)
    return(κms, κms_err)
end


#####################
## Compute the Q ratio structure to measure the clustering of the OC.
## See Cartwright et al. (2004)
## directly adapted to 3D
##

### !!!!!!! Normalization to be checked on λ , Cartwright & 2004
function get_Q(X, Y , Z)
    nxyz= length(X)
    A= Array{Float64}(undef,3,nxyz)
    for i in 1:nxyz
            A[:,i]= [X[i] Y[i] Z[i]]
    end

    d= Euclidean()
    p= pairwise(d, A, A, dims=2)

    src= Array{Int}(undef,0)
    dst= Array{Int}(undef,0)
    wgt= Array{Float64}(undef,0)

    for i in 1:nxyz
        for j in 1:nxyz
            if i != j
                push!(src,i)
                push!(dst,j)
                push!(wgt,p[i,j])
            end
        end
    end

    G= SimpleWeightedGraph(src, dst, wgt)
    kr= kruskal_mst(G, G.weights)

    ## compute the mean (arithmtic) edge of MST and of star separation. We do not apply other normalization
    λmst= 0
    for w in kr
        λmst+= w.weight
    end
    ## not the original normalization of Cartwright!
    k= (4π / 3 / nxyz )^(1/3)
    λmst /= k*length(kr)

    npair= nxyz*(nxyz-1)/2
    sbar= sum(p) / npair

    Q= λmst / sbar
    return(Q)
end
