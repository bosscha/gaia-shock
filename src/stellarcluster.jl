"""
functions and types to analyze stellar clusters using GAIA data

Author: s.leon @ ALMA

"""

### HR diagram metric using Delaunay tessalation to estimate density in the HRD 
###
## The tesselation algorithm works with coordinates in range [1,2]
    
function metricHRD(hrdi::Df , label)
    
    QA = []
    QP = []
    
    for ilab in label
        nhr = size(ilab)
        points = []
        
        for n in 1:nhr[1]
            push!(points,[hrdi.data[1,ilab[n]],hrdi.data[2,ilab[n]]])           
        end
    
        ## println("Voronoi tesselation ...")
        ## min_pts is the minimum of pts 
        if length(points) > 4
            per, area = voronoi(points)
        else
            per  = [1e10,1e10,1e10,1e10]
            area = [1e10,1e10,1e10,1e10]
        end
        
        p = sort(per ./ nhr[1])
        a = sort(area ./ nhr[1])
        
        n10 = max(1, trunc(Int, nhr[1] / 10))
        start10 = max(1 , 1*n10)
        end10   = min(nhr[1] , 9*n10)
        
        q1a  = mean(a[start10:end10])
        q1p  = mean(p[start10:end10])
        dq1a = std(a[start10:end10])
        dq1p = std(p[start10:end10])
        
        qa  = -log(max(1e-15, q1a))
        dqa = abs(dq1a/qa)   
        qp  = -log(max(1e-15, q1p)) 
        dqp = abs(dq1p / qp)
        
        push!(QA, (qa,dqa))
        push!(QP, (qp,dqp)) 
       
    end
    return(QA,QP)
end


##  proj :
##  spatial2d : y,z (id 2,3)
##  spatial3d : x,y,z (id 1,2,3)
##  velocity  : vrad,vdec (id 4,5)
##  HRD       : Diagram HR (id 6,7)
##
## labels : labels of the clusters##
##
## !!! works with dfcart !!! (cartesian and normalized)

function metric(s::Df, labels ,proj = "spatial2d", APERTURE = 1.0 , MAXAPERTURE = 15.0, NBOOTSTRAP = 50)

    if proj == "spatial2d"
        ycenter = mean(s.data[2,:])
        zcenter = mean(s.data[3,:]) 
    elseif proj == "spatial3d"
        xcenter = mean(s.data[1,:])      
        ycenter = mean(s.data[2,:])
        zcenter = mean(s.data[3,:]) 
    elseif proj == "velocity"
        vxcenter = mean(s.data[4,:])      
        vycenter = mean(s.data[5,:])
    elseif proj == "HRD"
        res = metricHRD(s, labels)
        return(res)
    end
    
    aper2 = APERTURE * APERTURE
    Q = []
    for ilab in labels
        slab = s.data[:,ilab]
            
        if proj == "spatial2d"
            yy = slab[2,:]
            zz = slab[3,:]
            ycenter = mean(yy)      
            zcenter = mean(zz)             
            angle_out = 2π * rand(Float64, NBOOTSTRAP)
            rad_out   = (MAXAPERTURE - 2APERTURE) * rand(Float64, NBOOTSTRAP)
            angle_in  = 2π * rand(Float64, NBOOTSTRAP)
            rad_in    = APERTURE * randexp(Float64, NBOOTSTRAP)
        elseif proj == "spatial3d"
            xx = slab[1,:]
            yy = slab[2,:]
            zz = slab[3,:]
            xcenter = mean(xx)             
            ycenter = mean(yy)      
            zcenter = mean(zz)             
            angle_out = 2π * rand(Float64, NBOOTSTRAP)
            phi_out   = π * rand(Float64, NBOOTSTRAP)
            rad_out   = (MAXAPERTURE - 2APERTURE) * rand(Float64, NBOOTSTRAP)
            angle_in  = 2π * rand(Float64, NBOOTSTRAP)
            phi_in    = π * rand(Float64, NBOOTSTRAP)           
            rad_in    = APERTURE * randexp(Float64, NBOOTSTRAP)
        elseif proj == "velocity"
            vx = slab[4,:]
            vy = slab[5,:]   
            vxcenter = mean(vx)      
            vycenter = mean(vy) 
            angle_out = 2π * rand(Float64, NBOOTSTRAP)
            rad_out   = (MAXAPERTURE - 2APERTURE) * rand(Float64, NBOOTSTRAP)
            angle_in  = 2π * rand(Float64, NBOOTSTRAP)
            rad_in    = APERTURE * rand(Float64, NBOOTSTRAP)   
        end  
            
        qc = zeros(NBOOTSTRAP)
            
        for k in 1:NBOOTSTRAP
            if proj == "spatial2d"
                yout = ycenter + rad_out[k] * cos(angle_out[k])
                zout = zcenter + rad_out[k] * sin(angle_out[k])
                radii_out = (slab[2,:] .- yout) .* (slab[2,:] .- yout) .+ (slab[3,:] .- zout) .* (slab[3,:] .- zout)
                yin = ycenter + rad_in[k] * cos(angle_in[k])
                zin = zcenter + rad_in[k] * sin(angle_in[k])
                radii_in = (slab[2,:] .- yin) .* (slab[2,:] .- yin) .+ (slab[3,:] .- zin) .* (slab[3,:] .- zin)                
            elseif proj == "spatial3d"
                xout = xcenter + rad_out[k] * cos(angle_out[k]) *  cos(phi_out[k])         
                yout = ycenter + rad_out[k] * sin(angle_out[k]) *  cos(phi_out[k])  
                zout = zcenter + rad_out[k] * sin(phi_out[k])
                radii_out = (slab[1,:] .- xout) .* (slab[1,:] .- xout) .+ (slab[2,:] .- yout) .* (slab[2,:] .- yout).+ (slab[3,:] .- zout) .* (slab[3,:] .- zout)
                xin = xcenter + rad_in[k] * cos(angle_in[k]) *  cos(phi_in[k])                  
                yin = ycenter + rad_in[k] * sin(angle_in[k]) *  cos(phi_in[k]) 
                zin = zcenter + rad_in[k] * sin(phi_in[k])
                radii_in = (slab[1,:] .- xin) .* (slab[1,:] .- xin) .+ (slab[2,:] .- yin) .* (slab[2,:] .- yin).+ (slab[3,:] .- zin) .* (slab[3,:] .- zin)              
            elseif proj == "velocity"
                vxout = vxcenter + rad_out[k] * cos(angle_out[k])
                vyout = vycenter + rad_out[k] * sin(angle_out[k])
                radii_out = (slab[4,:] .- vxout) .* (slab[4,:] .- vxout) .+ (slab[5,:] .- vyout) .* (slab[5,:] .- vyout)
                vxin = vxcenter + rad_in[k] * cos(angle_in[k])
                vyin = vycenter + rad_in[k] * sin(angle_in[k])
                radii_in = (slab[4,:] .- vxin) .* (slab[4,:] .- vxin) .+ (slab[5,:] .- vyin) .* (slab[5,:] .- vyin)                
            end
            
            isdout = radii_out .< aper2
            nout   = length(radii_out[isdout])
            isin   = radii_in .< aper2 
            nin    = length(radii_in[isin])
         
            qc[k] = log(max(1,nin)) - log(max(1,nout))
                
        end
            
        Qc    = mean(qc)
        Q_std = std(qc)    
        #println("Q: $Qc -- Q_std : $Q_std")
        push!(Q,(Qc,Q_std))  
    end
    return(Q)
end

### DBSCAN
function clusters(data , epsilon, leaf , minneigh, mincluster)
    eps = epsilon
    leafsize = leaf
    min_neighbors = minneigh
    min_cluster_size = mincluster

    res = dbscan(data , eps , leafsize = leaf, min_neighbors = minneigh, min_cluster_size=mincluster)
    
    label = Vector{Vector{Int}}()
    
    for cl in res
        indx = cl.core_indices
        append!(indx, cl.boundary_indices) 
        push!(label,indx)
    end
    return(label)
end