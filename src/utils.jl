### auxiliary methods
###

## return index which are not NAN
function  isnotnan(arr)
    index = map(!isnan, arr[:])
    
    return(index)
end

