## Main loop to extract in e2e the DBSCAN parameters.
## The list of votable is selected directly from the directory
##
## That version is to run in parallel the optimization

## OPTIONS
# -s : first centile to consider in the vot files (0-100)
# -e : last centile to consider in the vot files (0-100)
# -m : metafile for all split

using DataFrames, Query
using CSV, Glob, Dates
using Statistics, Random
using Printf


rootdir =  ENV["GAIA_ROOT"]

push!(LOAD_PATH,"$rootdir/run/src")
using GaiaClustering

## directory
wdir    = "$rootdir/e2e_products"
votdir  = "$wdir/votable"
plotdir = "$wdir/plotSelect"
ocdir   = "$wdir/oc"


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
    df         = filter_data(data, [0,2100])
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
        df= DataFrames.copy(CSV.read(blackname, delim=";"))
        blacklist= df.votname
    else
        blacklist= [""]
    end

    return(blacklist)
end

#######################
## Main loop
##
function main(filelist, metafile, prefile)
    let
        wd= pwd() ; nfile= size(filelist)[1] ; totalTime= 0.
        println("## Starting main loop using optimization with cycles...")
        println("## It can be very long but will be resumed to the last reduced file.")
        println("## Working directory: $wd")
        println("## $nfile files to analyze...")

        m= read_params(metafile, false)
        m.plotdir= plotdir
        m.ocdir= ocdir
        m.prefile= prefile
        fileres= "$(m.prefile).done.csv"

        # read a possible votname blacklist
        blackname= "blacklist.csv"
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
    print("## Main loop done.")
end
###############################################################################
let
    cd(votdir)
    votlist= glob("*.vot")
    cd(wdir)

    istart= 0 ; iend= 100
    metafile= "configAll.ext"

    for i in 1:length(ARGS)
        if ARGS[i] == "-s"
            istart= parse(Int, ARGS[i+1])
        end
        if ARGS[i] == "-e"
            iend= parse(Int, ARGS[i+1])
        end
        if ARGS[i] == "-m"
            metafile= ARGS[i+1]
        end
    end

    if istart < 0 istart= 0 end
    if iend > 100 iend= 100 end

###############
    header_extract()
    m= read_params(metafile, false)
    prefile= "$(m.prefile)-$istart-$iend"

    nfile= length(votlist)
    i1= convert(Int,floor(nfile*istart/100))
    if i1== 0 i1= 1 end

    if iend < 100
        i2= convert(Int,floor(nfile*iend/100)-1)
    else
        i2= convert(Int,floor(nfile*iend/100))
    end

    println("## DBSCAN optimization starting")
    println("## split files: $prefile")
    println("##")

    votsublist= votlist[i1:i2]
    main(votsublist, metafile, prefile)
    println("##\n## END of split optimization for $prefile \n##")
end
