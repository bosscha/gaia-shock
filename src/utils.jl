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
