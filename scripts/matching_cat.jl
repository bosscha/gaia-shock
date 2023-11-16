## Script to match two catalogues
## UGLY!!
using DataFrames , CSV , SQLdf

let 
    
TOLDIST= 0.1         ## distance tolerance (relative)
TOLANGLE= 0.05      ## angle tolerance on the sky in degrees

cat1= "CG2018_prepared.csv"
cat2= "catalogue_gaiadr3_final_bis.csv.merge"

# cat2= "spica_firda_milkyway_taufiq_alma_ocressc_joined_prepared.csv"

println("Matching $cat1 and $cat2 ...")
dat1 = CSV.read(cat1, DataFrame , delim=";")
dat2 = CSV.read(cat2, DataFrame , delim=";")

println(first(dat2))
new= true
nmatching= 0
df1= dat1 ; df2= dat2

for row1 in eachrow(dat1)
    dist1=row1.dmode
    ra1= row1.RAJ2000
    dec1= row1.DEJ2000
    nstar1= row1.Nstars

    for row2 in eachrow(dat2)
        dist2=row2.distance
        ra2= row2.ra
        dec2= row2.dec
        nstar2= row2.nstars

    
        if abs(ra1-ra2) < TOLANGLE && abs(dec1-dec2) < TOLANGLE && abs(dist1-dist2)/dist2 < TOLDIST
            dist2=row2.distance

            diff= 100*(dist1-dist2)/dist2
            ndiff= 100*(nstar1-nstar2)/nstar2
            println("diff: $diff - N: $nstar1 $nstar2")
            println("$ra1 $ra2 $dist1 $dist2")

            nmatching += 1
            if new 
                new= false
                df1 = similar(dat1,0)
                df2 = similar(dat2,0)
                push!(df1, row1)
                push!(df2, row2)
            else
                 push!(df1, row1)
                 push!(df2, row2)

            end
        end
    end
end

dmerge= hcat(df1,df2)

println("$nmatching matching...")
CSV.write("merge_cat.csv", dmerge, delim= ";")
end


