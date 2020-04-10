## check all votables and create a blacklist
##

using PyCall
using Distributions, Statistics , StatsBase

using Distances ,Random, Dates
using Glob

import PyPlot
sns= pyimport("seaborn")

## directory
rootdir = ENV["GAIA_ROOT"]
wdir    = "$rootdir/e2e_products"
votdir  = "$wdir/votable-TEST"
plotdir = "$wdir/plot-data-TEST"

push!(LOAD_PATH,"$rootdir/master/src")
using GaiaClustering

cd(wdir)

function plot_data(votname, df)

    s= size(df.data)
    npts= min(10000, s[2])
    idx= randperm(MersenneTwister(1234), s[2])[1:npts]

    PyPlot.plt.figure(figsize=(12.0,13.0))

    PyPlot.plt.subplot(2, 2, 1)
    PyPlot.plt.plot(df.data[1,idx],df.data[2,idx],"k.", markersize= 1)
    PyPlot.xlabel("l")
    PyPlot.ylabel("b")

    PyPlot.plt.subplot(2, 2, 2)
    PyPlot.plt.plot(df.raw[1,idx],df.raw[2,idx],"k.", markersize= 1)
    PyPlot.xlabel("α")
    PyPlot.ylabel("δ")

    PyPlot.plt.subplot(2, 2, 3)
    PyPlot.plt.plot(df.raw[1,idx],df.raw[5,idx],"k.", markersize= 1)
    PyPlot.xlabel("α")
    PyPlot.ylabel("parallax")

    PyPlot.plt.subplot(2, 2, 4)
    PyPlot.plt.plot(df.raw[6,idx],df.raw[7,idx],"k.", markersize= 1)
    PyPlot.xlabel("PM(α)")
    PyPlot.ylabel("PM(δ)")

    figname = plotdir*"/"*votname*".gaia.png"
    PyPlot.plt.savefig(figname)

    # PyPlot.plt.show()
end

function main(votlist)

    println("## Analyzing/plotting the raw data...\n")
    s=size(votlist)
    for i in 1:s[1]
        votname= votlist[i]
        data= read_votable(votdir*"/"*votname)
        df= filter_data(data)

        println("## Analyzing $votname")
        plot_data(votname,df)

    end
end

cd(votdir)
votlist= glob("*.vot")
cd(wdir)

main(votlist)
