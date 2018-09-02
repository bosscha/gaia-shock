## Run test with the Julia module gaiaClustering
## @02.09.2018

using gaiaClustering

## ... directories
rootdir = "/home/stephane/Science/cluster/GAIA"
wdir    = "/home/stephane/Science/cluster/GAIA/products"

cd(wdir)

voname = "NGC 2516-1.0deg.vot"

data       = read_votable(voname)
df         = filter_data(data)
dfcart     = add_cartesian(df)

blck       = [[1,2,3],[4,5], [6,7,8]]
wghtblck   = [3.,2.,1.5]
norm = "normal"
dfcartnorm = normalization_PerBlock(dfcart, blck, wghtblck , norm) 


println("## Test done...")
