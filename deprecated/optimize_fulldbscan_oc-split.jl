## Main loop to extract in e2e the DBSCAN parameters.
## The list of votable is selected directly from the directory
##
## That version is to run in parallel the optimization

## OPTIONS
# -s : first centile to consider in the vot files (0-100)
# -e : last centile to consider in the vot files (0-100)
# -o : root for the out file of the results

using DataFrames
using CSV, Glob, Dates
using Statistics
import DataFrames

rootdir =  ENV["GAIA_ROOT"]

push!(LOAD_PATH,"$rootdir/master/src")
using GaiaClustering

## directory
wdir    = "$rootdir/e2e_products"
votdir  = "$wdir/votable"
plotdir = "$wdir/plots-select-TEST"
ocdir   = "$wdir/oc-TEST"

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
    df         = filter_data(data)
    dfcart     = add_cartesian(df)
    blck       = [[1,2,3],[4,5], [6,7,8]]
    wghtblck   = [4.0,5.0,1.0]
    norm       = "identity"

    dfcartnorm , scale8 = normalization_PerBlock(dfcart, blck, wghtblck , norm, false)
    return(df, dfcart , dfcartnorm)
end

function mcmc_params()
    minQ    = 2.7
    minstars = 40
    forcedminstars = 30
##
    epsmean   = 2.5
    epsdisp   = 1.5
    min_nei   = 10
    min_cl    = 15
    ncoredisp = 10
    w3dmean   = 6.0
    w3ddisp   = 4.0
    wvelmean  = 6.0
    wveldisp  = 4.0
    whrdmean  = 2.0
    whrddisp  = 1.5
## MCMC parameters
    nburnout  = 2000
    niter     = 15000
##
    pinit = GaiaClustering.abcfull(minQ, minstars, forcedminstars, epsmean, epsdisp, min_nei, min_cl, ncoredisp, w3dmean, w3ddisp ,
    wvelmean, wveldisp, whrdmean, whrddisp, nburnout , niter)
     return(pinit)
end


function dbscanmcmcfull_updt!(ismcmc, fileres, mc ,votname)
    epsm = median(mc.eps)
    epsd = std(mc.eps)
    mneim = median(mc.mne)
    mneid = std(mc.mne)
    mclm = median(mc.mcl)
    mcld = std(mc.mcl)
    qcm = median(mc.qc)
    qnm = median(mc.qn)
    qcd = std(mc.qc)
    qnd = std(mc.qn)
    w3dm = median(mc.w3d)
    w3dd = std(mc.w3d)
    wvelm = median(mc.wvel)
    wveld = std(mc.wvel)
    whrdm = median(mc.whrd)
    whrdd = std(mc.whrd)

    println("## DBSCAN/MCMC stats:")
    println("### ϵ : ",epsm," +/- ", epsd)
    println("### min_nei  : ", mneim," +/- ", mneid)
    println("### min_clus  : ", mclm,"+/- ", mcld)
    println("### W3d  : ", w3dm,"+/- ", w3dd)
    println("### Wvel  : ", wvelm,"+/- ", wveld)
    println("### Whrd  : ", whrdm,"+/- ", whrdd)
    println("### Qn  : ",qnm," +/- ", qnd)
    println("### Qc  : ",qcm," +/- ", qcd)
    println("##")

    if !ismcmc
        res = DataFrame(votname=votname, epsm = epsm, epsd=epsd, mneim=mneim,mneid=mneid,mclm=mclm,mcld=mcld,
            qcm=qcm,qcd=qcd, qnm=qnm,qnd=qnd,
            w3dm=w3dm,w3dd=w3dd,wvelm=wvelm,wveld=wveld,whrdm=whrdm,whrdd=whrdd)
        CSV.write(fileres,res,delim=';')
        initmcmc = true
        return(res)
    else
        res = DataFrames.copy(CSV.read(fileres, delim=";"))
        newrow = DataFrame(votname=votname,epsm = epsm, epsd=epsd, mneim=mneim,mneid=mneid,mclm=mclm,mcld=mcld,
            qcm=qcm,qcd=qcd, qnm=qnm,qnd=qnd,
            w3dm=w3dm,w3dd=w3dd,wvelm=wvelm,wveld=wveld,whrdm=whrdm,whrdd=whrdd)
        println("### add DBSCAN/MCMC FULL results ...")

        append!(res,newrow)
        CSV.write(fileres,res,delim=';')
        return(newrow)
    end
end

function check_mcmc(votname, fileres)
    try
        res = DataFrames.copy(CSV.read(fileres, delim=";"))
        if votname in res.votname
            return(true, true)
        else
            return(false , true)
        end
    catch
        println("### $fileres will be created...")
        return(false, false)
    end
end

### update basic parameters of the extracted cluster
####
function SCparameters_updt(fileres,sc::GaiaClustering.SCproperties,votname)
    if isfile(fileres)
        res = DataFrames.copy(CSV.read(fileres, delim=";"))
        newrow = DataFrame(votname=votname)
        newrow = DataFrame(votname=votname,nstars=sc.nstars,distance=sc.distance,l=sc.l,b=sc.b,
            vl=sc.vl,vb=sc.vb,vrad=sc.vrad , xdisp=sc.xdisp,ydisp=sc.ydisp,zdisp=sc.zdisp,
            vldisp=sc.vldisp,vbdisp=sc.vbdisp , vraddisp=sc.vraddisp)
        println("### add SC basic properties...")
        append!(res,newrow)
        CSV.write(fileres,res,delim=';')
        return(newrow)
    else
        println("## No $fileres file, it will be created...")
        res = DataFrame(votname=votname,nstars=sc.nstars,distance=sc.distance,l=sc.l,b=sc.b,
            vl=sc.vl,vb=sc.vb,vrad=sc.vrad , xdisp=sc.xdisp,ydisp=sc.ydisp,zdisp=sc.zdisp,
            vldisp=sc.vldisp,vbdisp=sc.vbdisp , vraddisp=sc.vraddisp)
        CSV.write(fileres,res,delim=';')
        initmcmc = true
        return(res)
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
function main(filelist,fileres, fileSCres)
    let
        println("## Starting main loop ..")
        println("## It can be very long but will be resumed to the last reduced file.")
        wd= pwd()
        println("## Working directory: $wd")
        nfile= size(filelist)[1]
        println("## $nfile files to analyze...")
        totalTime= 0.

        # read a possible votname blacklist
        blackname= "blacklist-oc.csv"
        blacklist= read_blacklist(blackname)
        println("## blacklist read...")

        for i in 1:nfile
            votname = filelist[i]

            mcmcfound , ismcmcfile = check_mcmc(votname, fileres)

            ## test blacklist
            if votname in blacklist
                println("##")
                println("## $votname in blacklist...")
                println("##")
                mcmcfound= true
            end

            if !mcmcfound
                println("## Starting with $votname")
                tstart= now()
                println("## Starting at $tstart")

                df , dfcart , dfcartnorm = getdata(votdir*"/"*votname)

                ## MCMC optimization
                params = mcmc_params()
                mc = abc_mcmc_dbscan_full(dfcart, params)
                plot_dbscanfull_mcmc(plotdir, votname, mc , false)
                res = dbscanmcmcfull_updt!(ismcmcfile, fileres, mc ,votname)

                ## get the cluster and plot it
                println("## Extracting the cluster using DBSCAN/WEIGHTING with:")
                eps = res.epsm[1]
                min_nei = trunc(Int,res.mneim[1] + 0.5)
                min_cl = trunc(Int,res.mclm[1] + 0.5)
                w3d = res.w3dm[1]
                wvel = res.wvelm[1]
                whrd = res.whrdm[1]
                println("### ϵ : $eps")
                println("### min_neighbor: $min_nei")
                println("### min_cluster : $min_cl")
                println("### W3d : $w3d")
                println("### Wvel : $wvel")
                println("### Whrd : $whrd")

                mres = GaiaClustering.modelfull(eps,min_nei,min_cl,w3d,wvel,whrd)
                dfcartnorm = getDfcartnorm(dfcart, mres)
                labels = clusters(dfcartnorm.data ,eps  , 20, min_nei, min_cl)
                # labelmax , nmax = find_cluster_label(labels)
                labelmax , nmax = find_cluster_label2(labels, df, dfcart)
                println("### Label solution: $labelmax")
                println("### N stars: $nmax")

                ### Best solution...
                println("## Extracting and writing the OC stars...")
                export_df("$votname.0", ocdir, df , dfcart, labels , labelmax)
                scproperties0 = get_properties_SC(labels[labelmax] , df, dfcart)
                plot_cluster2(plotdir, "$votname.0", labels[labelmax], scproperties0,  dfcart , false)
                println("### ",scproperties0)

                #### all solutions are plotted/exported
                index=1
                for ilab in labels
                    scproperties = get_properties_SC(ilab , df, dfcart)
                    export_df("$votname.$index", ocdir, df , dfcart, labels , index)
                    plot_cluster2(plotdir, "$votname.$index", ilab, scproperties,  dfcart , false)
                    index +=1
                end
                ######

                SCparameters_updt(fileSCres, scproperties0, votname)

                tend= now()
                println("## Ending at $tend")
                duration= Dates.value(tend-tstart) / (1000*3600)
                totalTime += duration
                meanTime= totalTime / i
                ETA= meanTime * (nfile-i) / 24
                println("## Duration: $duration hours")
                println("## ETA: $ETA days")
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
    file_mcmc= "votlist-mcmc_full"

    for i in 1:length(ARGS)
        if ARGS[i] == "-s"
            istart= parse(Int, ARGS[i+1])
        end
        if ARGS[i] == "-e"
            iend= parse(Int, ARGS[i+1])
        end
        if ARGS[i] == "-o"
            file_mcmc= ARGS[i+1]
        end
    end

    if istart < 0 istart= 0 end
    if iend > 100 iend= 100 end

################
    filemcmcout= "$file_mcmc-$istart-$iend.csv"
    fileSCpropout= "$file_mcmc-SCproperties-$istart-$iend.csv"

    nfile= length(votlist)
    i1= convert(Int,floor(nfile*istart/100))
    if i1== 0 i1= 1 end

    if iend < 100
        i2= convert(Int,floor(nfile*iend/100)-1)
    else
        i2= convert(Int,floor(nfile*iend/100))
    end

    println("## DBSCAN optimization starting")
    println("## File out: $filemcmcout")
    println("## File SC properties: $fileSCpropout")

    votsublist= votlist[i1:i2]

    main(votsublist,filemcmcout ,fileSCpropout )

    println("##\n## END of split optimization...\n##")
end
