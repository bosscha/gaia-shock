## Testing an optimization on a small set of files
## It includes the optimization by cycle:detection+subtraction
##

## Main loop to extract in e2e the DBSCAN parameters.
## The list of votable is selected directly from the directory

using DataFrames, Query
using CSV, Glob, Dates
using Statistics, Random
using Printf

rootdir =  ENV["GAIA_ROOT"]

push!(LOAD_PATH,"$rootdir/run/src")
using GaiaClustering

## directory
wdir    = "$rootdir/products"
votdir  = "$rootdir/e2e_products/votable.2020"
plotdir = "$rootdir/products/test"
ocdir   = "$rootdir/products/octest"


## Maximum random votable for testing
MAX_VOTABLE = 10

## load a liist of votable and update the file if done
## add results
##

function readlist_votable(filelist::String)
    vot = DataFrames.copy(CSV.read(filelist))
    return(vot)
end

function getdata(filevot)
    voname = filevot

    data       = read_votable(voname)
    df         = filter_data(data, [0,2000])
    dfcart     = add_cartesian(df)
    blck       = [[1,2,3],[4,5], [6,7,8]]
    wghtblck   = [4.0,5.0,1.0]
    norm       = "identity"

    dfcartnorm , scale8 = normalization_PerBlock(dfcart, blck, wghtblck , norm, false)
    return(df, dfcart , dfcartnorm)
end


## to check if done and record
## check and updt if votname analyzed. If not done return false
function updt_votcompleted(fileres, votname , cycletot=1, flag= 0 , onlycheck=true)
    let
        if  onlycheck
            if !isfile(fileres)
                return(0, false)
            else
                res = DataFrames.copy(CSV.read(fileres, delim=";"))
                if votname in res.votname
                    x = @from i in res begin
                        @where i.votname == votname
                        @select i
                        @collect DataFrame
                    end
                    return(x, true)
                else
                    return(0, false)
                end
            end
            ## UPDTE
        else
            if !isfile(fileres)
                res = DataFrame(votname=votname, cycle=cycletot, flag=flag)
                CSV.write(fileres,res,delim=';')
                println("## $fileres created...")
                return(res,true)
            else
                res = DataFrames.copy(CSV.read(fileres, delim=";"))
                newrow = DataFrame(votname=votname,cycle=cycletot, flag=flag)
                append!(res,newrow)
                CSV.write(fileres,res,delim=';')
                return(res,true)
            end
        end
    end
end
## check for blacklist
##
function read_blacklist(blackname)
    if isfile(blackname)
        df= CSV.read(blackname, DataFrame, delim=";")
        blacklist= df.votname
    else
        blacklist= [""]
    end

    return(blacklist)
end

#######################
## Main loop
##
function main(filelist, metafile)
    let
        wd= pwd() ; nfile= size(filelist)[1] ; totalTime= 0.
        println("## Starting main loop using optimization with cycles...")
        println("## It can be very long but will be resumed to the last reduced file.")
        println("## Working directory: $wd")
        println("## $nfile files to analyze...")

        m= read_params(metafile, false)
        m.plotdir= plotdir
        m.ocdir= ocdir
        fileres= "$(m.prefile).done.csv"

        # read a possible votname blacklist
        blackname= "blacklist-test.csv"
        blacklist= read_blacklist(blackname)
        println("## blacklist read...")

        for i in 1:nfile
            votname = filelist[i]
            res, votfound= updt_votcompleted(fileres, votname, 0, 0, true)

            ## test blacklist
            if votname in blacklist
                println("##")
                println("## $votname in blacklist...")
                println("##")
                votfound= true
            end

            if !votfound
                tstart= now()
                println("###########################")
                println("## Starting with $votname")
                println("## Starting at $tstart")

                df , dfcart , dfcartnorm = getdata(votdir*"/"*votname)
                m.votname= votname
                cycle, flag= cycle_extraction(df, dfcart, m)
                res,  votfound= updt_votcompleted(fileres, votname , cycle, flag, false)

                tend= now()
                println("## Ending at $tend")
                println("## number of cycle: $cycle , flag:$flag ")
                println("##########################")
                println("##")

                duration= Dates.value(tend-tstart) / (1000*3600)
                totalTime += duration
                meanTime= totalTime / i
                ETA= meanTime * (nfile-i) / 24
                nleft= nfile-i
                ETAstr= @sprintf("%3.3f", ETA) ; durationstr= @sprintf("%3.3f", duration)
                @printf("## %s \n",specialstr("Duration: $durationstr hours","YELLOW"))
                @printf("## %s \n",specialstr("ETA: $ETAstr days","YELLOW"))
                @printf("## %s \n",specialstr("Votable done: $votname","YELLOW"))
                @printf("## %s \n",specialstr("Files analyzed: $i","YELLOW"))
                @printf("## %s \n",specialstr("Files to go: $nleft","YELLOW"))
                println("##\n##")

            end
        end
    end
    println("## Main loop done.")
end

###############################################################################
header_extract()

cd(votdir)
votlist= glob("*.vot")
cd(wdir)

rng = MersenneTwister()
shuffle!(rng, votlist)

main(votlist[1:MAX_VOTABLE],"configAll.ext")
