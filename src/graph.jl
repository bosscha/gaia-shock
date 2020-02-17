## Methods to be used with the Stellar Cluster associated with a Graph.
## Using Lightgraphs

## Zagreb-Randic indices
##

## first Zagreb index
## Normalized here with vertex nmber
function zagreb_first(g, α=2)
    nvertex= nv(g)
    deg= degree(g)
    M1= sum(deg .^ α) / nvertex

    return(M1)
end


## second Zagreb index
##
function zagreb_second(g, α=1 )
    nvertex= nv(g)
    nedge= ne(g)
    edg= edges(g)

    M2= 0
    for e in edg
        d1= degree(g,src(e))
        d2= degree(g,dst(e))
        M2 += (d1*d2)^α
    end

    M2 /= nedge

    return(M2)
end
