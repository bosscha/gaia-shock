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

     dct= Dict()
     for field in f
         val= getfield(s, field)
         push!(dct, field => val)
     end

     df= DataFrame(dct)
     return(df)
end
