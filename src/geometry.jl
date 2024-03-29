### Function to analyze the geometry
###
function voronoi_perimeter(ver, region)
    # compute the perimetet  in the voronoi region and vertices

    perimeter = 0.
    sizeRegion = length(region)

    if -1 in region
        perimeter = 0.
    else
        for index in 1:sizeRegion
            index2 = (index%sizeRegion) + 1

            if region[index]!= -1 || region[index2] != -1

                x0 = ver[region[index]+1,1]
                y0 = ver[region[index]+1,2]

                x1 = ver[region[index2]+1,1]
                y1 = ver[region[index2]+1,2]
                perimeter += sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
            end
        end
    end
    return(perimeter)
end
##
function voronoi_area(ver, region)
    # compute the area   in the voronoi region and vertices

    area = 0.
    sizeRegion = length(region)

    if -1 in region
        area = 0.
    else
        for index in 1:sizeRegion
            x0 = ver[region[index]+1,1]
            y0 = ver[region[index]+1,2]

            index2 = (index%sizeRegion) + 1
            x1 = ver[region[index2]+1,1]
            y1 = ver[region[index2]+1,2]

            area += 0.5 * (x0*y1 - x1*y0)
        end
    end
    return(abs(area))
end
## Voronoi tesselation using the scipy function
##
function voronoi_python(pts, verbose = true)
    let
        ndat = length(pts)
        peri = zeros(ndat)
        area = zeros(ndat)

        vor = 0
        spatial = pyimport("scipy.spatial")

        try
            # vor = spatial[:Voronoi](pts)
            vor = spatial.Voronoi(pts)
        catch
            println("## Voronoi error...")
            return([1000,1000],[1000,1000])   ## arbitrary values for peri and area
        end

        if verbose println("## Voronoi tesselation done.") end

        ver = vor.vertices
        reg = vor.regions
        pt  = vor.point_region

        for i in 1:ndat
            region  =  reg[pt[i]+1]
            peri[i] =  voronoi_perimeter(ver,region)
            area[i] =  voronoi_area(ver,region)
        end
        println("ending voronoi ...")
        return(peri , area)
    end
end
##
## Julia voronoi implementation
function voronoi(pts, verbose = false)

        ndat = length(pts)
        peri = zeros(ndat)
        area = zeros(ndat)

        xx= [] ; yy=[]
        for i in 1:ndat
            push!(xx,pts[i][1])
            push!(yy,pts[i][2])
        end

        xmin= minimum(xx)
        xmax= maximum(xx)
        ymin= minimum(yy)
        ymax= maximum(yy)

        rect = Rectangle(gb.Point2(xmin,ymin),gb.Point2(xmax,ymax))
        ptarr= [gb.Point2(xx[i], yy[i]) for i in 1:ndat]

        tess = voronoicells(ptarr, rect)
        area= voronoiarea(tess)

        return(area , area)
end
