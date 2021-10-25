## info on the current run for the extraction
## analysis of the .csv produced

using DataFrames, Query
using CSV, Glob, Dates
using Printf

function perf_status(file)
    if isfile(file)
        perf= CSV.File(file, delim=";")
        nstarcluster= sum(perf.nmax)

        x = @from i in perf begin
            @where i.cycle == 1
            @select i
            @collect DataFrame
        end
        nstartotal= sum(x.nstar)
    else
        nstarcluster= 0
        nstartotal= 0
    end
    return(nstarcluster, nstartotal)
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

function bold(s) specialstr(s,"BOLD") end
function yellow(s) specialstr(s,"YELLOW") end

function main(csvlist)
    println("## Extraction status")
    println("## Based on files *.done.csv and *.perf.csv")
    println("## $(now())")
    println("##")

    nvotTotal= 0 ; ncycleTotal= 0 ; nstarclusterTotal= 0 ; nstarTotal= 0
    for file in csvlist
        done= CSV.File(file, delim=";")

        nvot= length(done.votname)
        ncycle= sum(done.cycle)
        fileperf= split(file,".")[1]*".perf.csv"
        nstarcluster, nstartotal= perf_status(fileperf)

        println("### $file")
        println("### votables: $nvot")
        println("### cycles: $ncycle")
        println("### stars extracted: $nstarcluster")
        println("### total star analyzed: $nstartotal")
        println("###")

        nvotTotal += nvot
        ncycleTotal += ncycle
        nstarclusterTotal += nstarcluster
        nstarTotal += nstartotal
    end
    println("##")
    println(yellow("## SUMMARY:"))
    println(yellow("## votables: $nvotTotal"))
    println(yellow("## cycles: $ncycleTotal"))
    println(yellow("## stars extracted: $nstarclusterTotal"))
    println(yellow("## stars analyzed: $nstarTotal"))

end

list= glob("*.done.csv")
main(list)
