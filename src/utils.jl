### auxiliary methods
###

## return index which are not NAN
function  isnotnan(arr)
    index = map(!isnan, arr[:])

    return(index)
end

###return blacklist if any
function read_blacklist(blackname)
    if isfile(blackname)
        df= CSV.read(blackname, delim=";")
        blacklist= df.votname
    else
        blacklist= [""]
    end

    return(blacklist)
end

## transform struct to DF
function convertStruct2Df(s)::DataFrame
    T= typeof(s)
    f= fieldnames(T)
    n=length(f)

     field1= fieldname(T,1) ; val1= getfield(s, field1)
     df= DataFrame(field1 => val1)

     for i in 2:n
         field= fieldname(T,i)
         val= getfield(s, field)
         insertcols!(df,i,field=>val)
     end

     return(df)
end

## special printing code
## return string
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
function purple(s) GaiaClustering.specialstr(s,"PURPLE") end
function cyan(s) GaiaClustering.specialstr(s,"CYAN") end
function blue(s) GaiaClustering.specialstr(s,"BLUE") end
function red(s) GaiaClustering.specialstr(s,"RED") end

function header_extract()
    println("#")
    println("# >>>>>>>>>> EXTRActing Stellar Clusters from GAIA data <<<<<<<<<<")
    println("# Version: $VERSION  (https://github.com/bosscha/gaia-shock)")
    println("#\n#")
end

function debug_red(msg)
    # println(red("##_debug_###### $msg"))
end
