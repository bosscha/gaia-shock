## extract the final open clusters selection with the DBSCAN parameters
##

using PyCall
using Distributions
using Statistics
using DataFrames

using Distances
using Random , Printf

rootdir = ENV["GAIA_ROOT"]

push!(LOAD_PATH,"$rootdir/master/src")
using GaiaClustering

import PyPlot , CSV

wdir    = "$rootdir/e2e_products"
votdir  = "$rootdir/e2e_products/votable"
ocdir   = "$rootdir/e2e_products/oc"
cd(wdir)

fileparams= "votlist.finalSample_metric0.020.csv"

###########################################################################
function getdata(filevot, wghtblck)
    voname     = filevot
    data       = read_votable(voname)
    df         = filter_data(data)
    dfcart     = add_cartesian(df)
    blck       = [[1,2,3],[4,5], [6,7,8]]
    norm       = "identity"

    dfcartnorm , scale8 = normalization_PerBlock(dfcart, blck, wghtblck , norm, false)
    return(df, dfcart , dfcartnorm)
end

## Check if OC was extracted
function _check_extraction(votname, fileres)
    try
        res = CSV.read(fileres, delim=";")
        if votname in res.votname
            return(true)
        else
            return(false)
        end
    catch
        println("### $fileres will be created...")
        return(false)
    end
end

function _extraction_updt(fileres, votname)
    try
        res = CSV.read(fileres, delim=";")
        newrow = DataFrame(votname=votname)
        append!(res,newrow)
        CSV.write(fileres,res,delim=';')
    catch
        println("## No $fileres file, it will be created...")
        res= DataFrame(votname=votname)
        CSV.write(fileres,res,delim=';')
    end
end

## Create the DataFrame to save the cluster...
##
function _export_df(votname, ocdir, df , dfcart, labels , labelmax)
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

    oc= DataFrame(ra=ra,dec=dec,l=l,b=b, distance=d,pmra=pmra, pmdec=pmdec, X=X,Y=Y,Z=Z,vl=vl,vb=vb,vrad=vrad,gbar=gbar,rp=rp,bp=bp, ag=ag)

    filename= @sprintf("%s/%s-oc.csv",ocdir, votname[1:end-4])
    CSV.write(filename,oc,delim=';')
    @printf("\n## %s created",filename)
end

## Main loop
##

function main(paramfile, fileres)
    let
        println("## Starting main loop..")
        spl= CSV.read(paramfile, delim= ";")
        println("## $paramfile read..")
        println("## Check if csv file is ; separated..")
        s=size(spl)

        for i in 1:s[1]
            votname = spl.votname[i]
            println("## Starting with $votname")
            wght= [spl.w3dm[i],spl.wvelm[i] ,spl.whrdm[i]]
            println(wght)

            found= _check_extraction(votname, fileres)
            if !found
                df , dfcart , dfcartnorm = getdata(votdir*"/"*votname, wght)

                ## get the cluster and plot it
                println("## Extracting the cluster using DBSCAN with:")
                eps = spl[:epsm][i]
                min_nei = trunc(Int,spl[:mneim][i] + 0.5)
                min_cl = trunc(Int,spl[:mclm][i] + 0.5)
                println("### ϵ : $eps")
                println("### min_neighbor: $min_nei")
                println("### min_cluster : $min_cl")
                labels = clusters(dfcartnorm.data ,eps  , 20, min_nei, min_cl)
                labelmax , nmax = find_cluster_label(labels)
                println("### Label solution: $labelmax")
                println("### N stars: $nmax")

                _export_df(votname, ocdir, df , dfcart, labels , labelmax)
                _extraction_updt(fileres, votname)
            end
        end
    end
end

main(fileparams , "votlist.finalSample.extracted.csv")
