## install package for Gaia.Clustering

using Pkg
ENV["PYTHON"] = " "

Pkg.add("PyCall")
Pkg.build("PyCall")

Pkg.add("DataFrames")
Pkg.add("Clustering")
Pkg.add("Statistics")
Pkg.add("Distributions")
Pkg.add("Random")
Pkg.add("CSV")
Pkg.add("PyPlot")
Pkg.add("Formatting")
Pkg.add("Query")
Pkg.add("TSne")
Pkg.add("Printf")
Pkg.add("Interpolations")
Pkg.add("Images")
Pkg.add("LightGraphs")
Pkg.add("SimpleWeightedGraphs")
Pkg.add("Distances")
Pkg.add("Glob")
Pkg.add("StatsBase")
Pkg.add("FFTW")
Pkg.add("LaTeXStrings")
Pkg.add("MultivariateStats")
Pkg.add("ArgParse")
Pkg.add("VoronoiCells")
Pkg.add("GeometryBasics")

println("## Package installation for GaiaClustering done...")
