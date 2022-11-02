## functions to perform imaging computation
##

b3spline1d = [1/16 , 1/4 , 3/8 , 1/4 ,  1/16]

## Wavelet A Trous Algorithm (See e.g. Leon et al. 2001)
## im: image array
## levels : number of Wavelet planes
## kernel1d : 1D kernel to create the 2D kernel. Default is B3-spline but  could be other (e.g. [0,1,0])
##
function atrous(im, levels , kernel1d = b3spline1d, verbose = true)
    let
        w      = []
        ker    = Base.copy(kernel1d)
        approx = Base.copy(im)
         
        for lev in 1:levels    
            if verbose
                println("## WT--A Trous--Plane: $lev")
            end    
        
            kernel2d = zeros(length(ker),length(ker))
            for i in 1:(length(ker))
                kernel2d[i,:] = ker[i] .* ker
            end
            ker2 = ImageFiltering.centered(kernel2d)
       
            smooth = ImageFiltering.imfilter(approx, ker2)  
            wi = approx .- smooth
            approx = smooth
        
            ker = upscale(ker)
            push!(w, wi) 
        end
    
        ## LSP ###
        if verbose
            println("## WT--A Trous--LSP ($(levels+1))")
        end
        
        push!(w, approx)
    
        return(w)
    end
end

## upscale by a factor 2 the 1D kernel 
function upscale(kerr)
    xx   = 0:length(kerr)-1
    itp  = Interpolations.interpolate(kerr, Interpolations.BSpline(Interpolations.Linear()))
    sitp = Interpolations.scale(itp, xx)
    xp   = 0:0.5:length(kerr)-1
    ker2 = 0.5 .* sitp(xp)
    
    return(ker2) 
end

### Add the WT coefficients to get the (partial) image
##
function addWav(wt,p1,p2)
    let
        
        nplane = size(wt)
        if p2 == -1 p2 = nplane[1] end
        nxy = size(wt[1])
        println(nxy)
        im = zeros(nxy[1], nxy[2])
    
        for i in p1:p2
            im = im .+ wt[i]
        end
    return(im)
    end
end

## Hard thresholding on the WT coefficients.
## wt: wavelet transform
## d : random distribution
## σ : threshold in std from the simulated noise, if an array of length== nw[1] then used for thresholding
##
function thresholdingWav(wt, d::UnivariateDistribution , σ = 3 , kernel1d=b3spline1d)
    wtfilt = Base.copy(wt)
    nw = size(wt)
    nxy = size(wt[1])
    
    println("## WT--Thresholding the WT")
    if length(σ) == nw[1]
        threshold = σ
    else
        println("## WT--Simulating the noise in the WT.")
        threshold = zeros(nw[1])
        noise = rand(d , nxy[1] , nxy[2])
        wsimul = atrous(noise , nw[1]-1 , kernel1d , false)
        threshold = σ .* noiseWav(wsimul)
    end
    
    for p in 1:nw[1]
        for i in 1:nxy[1]
            for j in 1:nxy[2]
                if abs(wt[p][i,j]) < threshold[p]
                   wtfilt[p][i,j] = 0
                end
            end
        end
    end
    return(wtfilt)
end


## Compute std per plane for the wt
##
function noiseWav(wt)
    noise = []
    for plane in wt
        σ = std(plane)
        push!(noise,σ)
    end
    return(noise) 
end

