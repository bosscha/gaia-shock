## Filter the MCMC results for DBSCAN on the ϵ+Xdisp metric (see NBs)
##

using DataFrames , CSV , Statistics
using Query
using PyPlot
using TSne
using Printf

import CSV

rootdir = ENV["GAIA_ROOT"]
wdir    = "$rootdir/e2e_products"
plotdir = "$rootdir/e2e_products"

cd(wdir)

## Input dataframe
mcfull= CSV.read("votlist-mcmc_full-TEST.csv", delim = ";")
scfull= CSV.read("votlist-SCproperties_full-TEST.csv", delim = ";")
dfcomb= join(mcfull, scfull, on = :votname, makeunique=true)

## Building the metric
x = @from i in dfcomb begin
    @let votname=i.votname
    @let q1 = (i.epsm / i.mclm) * (i.xdisp / i.distance)

    @select {votname=votname,Q1=q1}
    @collect DataFrame
end

println(x)

## Plot the metric histogram
fig = figure(figsize=(7,5))
nbins=50
h = plt.hist(x.Q1,nbins,range = [0,0.15], color = "g",  label = "dbscan+weight density" , alpha=0.6 , density =true)
plt.xlabel("ϵ+Xdisp Metric")
plt.savefig("dbscanMetric.png")

## filtering
# Threshold maximum on Q1 to filter the OC. Check it in the plot dbscanMetric
#
THRESHOLD_Q1= 0.02

final = @from i in dfcomb begin
    @where (i.xdisp / i.distance)*(i.epsm / i.mclm) < THRESHOLD_Q1
    @select i
    @collect DataFrame
end

println("# N selected")
println(size(final))
fileout= @sprintf("votlist.finalSample_metric%3.3f.csv", THRESHOLD_Q1)
CSV.write(fileout, final, delim = ';')
