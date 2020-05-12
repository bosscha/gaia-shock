## functions to set/read the different parameters for the extraction
##
## return the default parameters
##
function set_default_params()::meta
    def = meta()

    ## Directories
    def.rootdir =  ENV["GAIA_ROOT"]
    def.wdir    = @sprintf("%s/e2e_products", def.rootdir)
    def.votdir  = @sprintf("%s/votable", def.wdir)
    def.plotdir = @sprintf("%s/plots-select" , def.wdir)
    def.ocdir   = @sprintf("%s/oc", def.wdir)
    def.votname = "test.vot"

    ## MCMC
    ##
    def.minQc    = 2.7
    def.minQn    = 40
    def.maxQn    = 10000
    def.forcedminstars = 30
    def.mingoodsolution = 10        # minimum good solution for Q in check_qminqstar
    def.niterqminq = 500            # check_qminqstar..

    def.nburnout  = 2000            # burn-in iterations
    def.nchain    = 15000           # MCMC iterations==number of CHAINS
    def.maxiter   = 500000


    ## prior settings
    ##
    def.epsmean   = 2.5             # epsilon mean
    def.epsdisp   = 1.5             # epsilon dispersion
    def.min_nei   = 10              # min_nei mean
    def.min_cl    = 15              # min_cl mean
    def.ncoredisp = 10              # dispersion for min_nei and min_cl
    def.w3dmean   = 6.0
    def.w3ddisp   = 4.0
    def.wvelmean  = 6.0
    def.wveldisp  = 4.0
    def.whrdmean  = 2.0
    def.whrddisp  = 1.5

    ## Metrics find_clusters
    def.aperture2d    = 1.5
    def.maxaperture2d = 15
    def.aperture3d    = 3
    def.maxaperture3d = 20
    def.aperturev     = 3.0
    def.maxaperturev  = 20
    def.nboot         = 50          # bootstrap number

    def.labels        = "Qc"        # method to get labels: Qc|Qn

    return(def)
end

## extracting meta from .ext file
function read_params(file)
    par= set_default_params()
    println("## All parameters set to default...")
    t= readdlm(file,comments=true, comment_char='#')

    s= size(t)
    for i in 1:s[1]
        println(t[i,:])
        if t[i,2] == "="
            set_param!(par,t[i,1],t[i,3])
        else
            p= t[i,1]
            println("### parsing parameter file could be wrong near $p ...")
        end
    end

    println("## Parameters read from $file")

    return(par)
end

## for parsing .ext file
function set_param!(def, parstr,value)
    if parstr == "rootdir" def.rootdir= value end
    if parstr == "wdir"    def.wdir= value end
    if parstr == "votdir"  def.votdir= value end
    if parstr == "plotdir" def.plotdir= value end
    if parstr == "ocdir"   def.ocdir= value end
    if parstr == "votname" def.votname= value end

    if parstr == "minQc"           def.minQc= value end
    if parstr == "minQn"           def.minQn= value end
    if parstr == "maxQn"           def.maxQn= value end
    if parstr == "forcedminstars"  def.forcedminstars= value end
    if parstr == "mingoodsolution" def.mingoodsolution= value end
    if parstr == "niterqminq"      def.niterqminq= value end

    if parstr == "nburnout" def.nburnout= value end
    if parstr == "nchain"   def.nchain= value end
    if parstr == "maxiter"  def.maxiter= value end

    if parstr == "epsmean"   def.epsmean= value end
    if parstr == "epsdisp"   def.epsdisp= value end
    if parstr == "min_nei"   def.min_cl= value end
    if parstr == "min_cl"    def.min_cl= value end
    if parstr == "ncoredisp" def.ncoredisp= value end
    if parstr == "w3dmean"   def.w3dmean= value end
    if parstr == "w3ddisp"   def.w3ddisp= value end
    if parstr == "wvelmean"  def.wvelmean= value end
    if parstr == "wveldisp"  def.wveldisp= value end
    if parstr == "whrdmean"  def.whrdmean= value end
    if parstr == "whrddisp"  def.whrddisp= value end

    ## Metrics find_clusters
    if parstr == "aperture2d"     def.aperture2d= value end
    if parstr == "maxaperture2d"  def.maxaperture2d= value end
    if parstr == "aperture3d"     def.aperture3d= value end
    if parstr == "maxaperture3d"  def.maxaperture3d= value end
    if parstr == "aperturev"      def.aperturev= value end
    if parstr == "maxaperturev"   def.maxaperturev= value end
    if parstr == "nboot"          def.nboot= value end

    if parstr == "labels"  def.labels= value end

end
