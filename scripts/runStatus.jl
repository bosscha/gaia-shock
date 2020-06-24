## info on the current run for the extraction
## analysis of the .csv produced

using DataFrames, Query
using CSV, Glob, Dates

function perf_status(file)
    if isfile(file)
        perf= DataFrames.copy(CSV.read(file, delim=";"))
        nstarcluster= sum(perf.nmax)
    else
        nstarcluster= 0
    end
    return(nstarcluster)
end

function specialstr(str, CODE)
    d= Dict( "PURPLE" => "\033[95m" ,
             "CYAN" => "\033[96m" ,
             "DARKCYAN" => "\033[36m" ,
             "BLUE" => "\033[94m" ,
             "GREEN" => "\033[92m" ,
             "YELLOW" => "\033[93m" ,
             "RED" => "\033[91m" ,
             "BOLD" => "\033[1m" ,
             "UNDERLINE" => "\033[4m" ,
             "END" => "\033[0m")

    if haskey(d, CODE)
        res= @sprintf("%s%s%s", d[CODE], str, d["END"])
        return(res)
    else
        return(str)
    end
end

function bold(s) GaiaClustering.specialstr(s,"BOLD") end
function yellow(s) GaiaClustering.specialstr(s,"YELLOW") end

function main(csvlist)
    println("## Extraction status")
    println("## Based on files *.done.csv and *.perf.csv")
    println("## $(now())")
    println("##")

    nvotTotal= 0 ; ncycleTotal= 0 ; nstarclusterTotal= 0
    for file in csvlist
        done= DataFrames.copy(CSV.read(file, delim=";"))

        nvot= length(done.votname)
        ncycle= sum(done.cycle)
        fileperf= split(file,".")[1]*".perf.csv"
        nstarcluster= perf_status(fileperf)

        println("### $file")
        println("### votables: $nvot")
        println("### cycles: $ncycle")
        println("### stars extracted: $nstarcluster")
        println("###")

        nvotTotal += nvot
        ncycleTotal += ncycle
        nstarclusterTotal += nstarcluster
    end
    println("##")
    println(yellow("## SUMMARY:"))
    println(yellow("## votables: $nvotTotal"))
    println(yellow("## cycles: $ncycleTotal"))
    println(yellow("## stars extracted: $nstarclusterTotal"))

end

list= glob("*.done.csv")
main(list)
